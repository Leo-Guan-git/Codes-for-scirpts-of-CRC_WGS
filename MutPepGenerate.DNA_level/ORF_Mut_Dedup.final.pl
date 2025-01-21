#! /usr/bin/env perl -w
# @Author: Xiangyu Guan
# @Date:   2025-01-21
# @Last Modified by:   Xiangyu Guan
# @Last Modified time: 2025-01-21
use strict;
use warnings;
use POSIX;
use File::Basename;
use Data::Dumper;

@ARGV>=4 or die "Usage: perl $0 <ference fasta file> <output tsv file> <output fasta file> <input file[s]>\n";
my ($ref_f, $out_tsv_f, $out_fa_f, @inputs) = @ARGV;
my (%tmp, %contents, %remained, %strings, %origin, %Refs);
my ($headline, $max_len, $min_len, $check, $out_tmp_prefix);
my %AA_condon = ('GCT'=>'A', 'GCC'=>'A', 'GCA'=>'A', 'GCG'=>'A', 'CGT'=>'R', 'CGC'=>'R',
                 'CGA'=>'R', 'CGG'=>'R', 'AGA'=>'R', 'AGG'=>'R', 'TCT'=>'S', 'TCC'=>'S',
                 'TCA'=>'S', 'TCG'=>'S', 'AGT'=>'S', 'AGC'=>'S', 'ATT'=>'I', 'ATC'=>'I',
                 'ATA'=>'I', 'TTA'=>'L', 'TTG'=>'L', 'CTT'=>'L', 'CTC'=>'L', 'CTA'=>'L',
                 'CTG'=>'L', 'GGT'=>'G', 'GGC'=>'G', 'GGA'=>'G', 'GGG'=>'G', 'GTT'=>'V',
                 'GTC'=>'V', 'GTA'=>'V', 'GTG'=>'V', 'ACT'=>'T', 'ACC'=>'T', 'ACA'=>'T',
                 'ACG'=>'T', 'CCT'=>'P', 'CCC'=>'P', 'CCA'=>'P', 'CCG'=>'P', 'AAT'=>'N',
                 'AAC'=>'N', 'GAT'=>'D', 'GAC'=>'D', 'TGT'=>'C', 'TGC'=>'C', 'CAA'=>'Q',
                 'CAG'=>'Q', 'GAA'=>'E', 'GAG'=>'E', 'CAT'=>'H', 'CAC'=>'H', 'AAA'=>'K',
                 'AAG'=>'K', 'TTT'=>'F', 'TTC'=>'F', 'TAT'=>'Y', 'TAC'=>'Y', 'ATG'=>'M',
                 'TGG'=>'W', 'TAG'=>'*', 'TGA'=>'*', 'TAA'=>'*'
                 );
my $output_dir = dirname $out_tsv_f;
sub Message{
    my $message = $_[0];
    my $datestring = localtime();
    print "\t$datestring >>> $message\n";
}

Message("DeDup Job Start!");
$out_tmp_prefix = $output_dir."/tmp";
if (! -e $output_dir){ system("mkdir -p $output_dir"); }

foreach my $in_f(@inputs){
    $in_f =~ /\.PEP(\d+)\./;
    $origin{$1} = $in_f;
}

my @length = (sort {$b <=> $a} keys %origin);
$max_len = $length[0];
$min_len = $length[$#length];
Message("Max Length of Peptides: $max_len");
Message("Min Length of Peptides: $min_len");
my $sub_title;
foreach my $len(sort {$b <=> $a} @length){
    my $in_f = $origin{$len};
    $sub_title += 1;
    # Message("$in_f Dedup start...");
    # %tmp = ();
    # %strings = ();
    %Refs = ();
    Message("Step 1.${sub_title}: reference database $ref_f $len-mer Parsing");
    # %Refs = RefProIndexing($ref_f, $len);
    RefProIndexing($ref_f, $len);
    Message("Step 1.${sub_title} Done!");
    Message("Step 2.${sub_title}: input file $in_f remove Redundance");
    foreach my $ref_len(sort {$b <=> $a} $len+1..$max_len){
        my $tmp_ref_f = "$out_tmp_prefix.PEP${ref_len}.final.tsv";
        Message("Step 2.${sub_title}: $tmp_ref_f $len-mer Parsing");
        if( ! -e $tmp_ref_f ){ next; }
        open REF, "< $tmp_ref_f";
        readline REF;
        while(my $line=<REF>){
            chomp $line;
            my $pep = (split /\t/, $line)[9];
            my $kmer_num = $ref_len - $len + 1;
            for (my $n = 0; $n < $kmer_num; $n++) {
                my $ind = substr $pep, $n, $len;
                $Refs{$ind} = 1;
            }
        }
        close REF;
    }
    my $out_tmp_f = "$out_tmp_prefix.PEP$len.final.tsv";
    if($in_f =~ /\.gz$/){ open IN, "gzip -dc $in_f |" or die $!; }else{ open IN, "< $in_f" or die $!; }
    # $headline= Dedup($in_f);
    $headline = readline IN;
    %remained = ();
    while (my $line=<IN>){
        chomp $line;
        my $synon = Synonymous($line);
        if($synon == 1){ next; }
        my $PEP = (split /\t/, $line)[9];
        if($Refs{$PEP} == 1){ next; }
        unless (exists $remained{$PEP}){
            $remained{$PEP} = $line;
        }
    }
    close IN;
    Writetmp($headline, $out_tmp_f);
    Message("Step 2.${sub_title} Done!");
}

Message("step3: Result write to $out_tsv_f and $out_fa_f");

if ($out_tsv_f =~ /\.gz$/){ open TSV, "| gzip -c > $out_tsv_f " or die $!; }else{ open TSV, "> $out_tsv_f" or die $!; }
if($out_fa_f =~ /\.gz$/){ open FA, "| gzip -c > $out_fa_f" or die $!; }else{ open FA, "> $out_fa_f" or die $!; }

if($ref_f =~ /\.gz$/){ open REF, "gzip -dc $ref_f |" or die $!; }else{ open "< $ref_f" or die $!; }
my $linebytes = 0;
while(my $line=<REF>){
    print FA $line;
    chomp $line;
    $linebytes = length($line) if(($line !~ /^>/) && (length($line) > $linebytes));
}
close REF;

my ($MP_num, $pseudo_num, $pseudo_len, $pseudo_seq) = (0,0,500,"");
my $break = 'B' x 10;
my @Included;
print TSV "MP_ID\t$headline";
foreach my $len(sort {$b <=> $a} @length){
    my $tmp_in_f = "$out_tmp_prefix.PEP$len.final.tsv";
    open IN, "< $tmp_in_f";
    readline IN;
    while(<IN>){
        chomp;
        print TSV "MP${MP_num}\t$_\n";
        push @Included, $MP_num;
        $MP_num +=1;
        my $PEP = (split /\t/, $_)[9];
        $pseudo_seq .= $PEP.$break;
        if(length($pseudo_seq) >= $pseudo_len){
            my ($begin, $end) = (sort {$a <=> $b} @Included)[0,$#Included];
            print FA ">MUTPRO|MUTPRO${pseudo_num}|MUTPRO${pseudo_num} MP${begin}-${end}\n";
            for (my $n = 0; $n < length($pseudo_seq); $n+=$linebytes) {
                print FA substr($pseudo_seq, $n, $linebytes);
                print FA "\n";
            }
            @Included = ();
            $pseudo_num +=1;
            $pseudo_seq = "";
        }
    }
    close IN;
    # system("rm $tmp_in_f") if (-e "$tmp_in_f");
}
if($pseudo_seq){
    my ($begin, $end) = (sort {$a <=> $b} @Included)[0,$#Included];
    print FA ">MUTPRO|MUTPRO${pseudo_num}|MUTPRO${pseudo_num} MP${begin}-${end}\n";
    for (my $n = 0; $n < length($pseudo_seq); $n+=$linebytes) {
        print FA substr($pseudo_seq, $n, $linebytes)."\n";
    }
    $pseudo_num +=1;
}
close TSV;
close FA;
Message("step3 Done");
Message("Selected $MP_num mutated peptides");
Message("Generated $pseudo_num pseudo-proteins");

sub NS2PEP{
    my ($ns_seq, $strand, $aa_len) = @_;
    my $pep_seq;
    if ($strand eq "-"){
        $ns_seq =~ tr/ATCG/TAGC/;
        $ns_seq = reverse $ns_seq;
    }
    if(length($ns_seq) < $aa_len*3){ return ""; }
    $ns_seq = substr $ns_seq, 0, $aa_len*3;
    for (my $n = 0; $n < length($ns_seq); $n += 3) {
        my $condon = substr $ns_seq, $n, 3;
        $pep_seq .= $AA_condon{$condon};
    }
    return $pep_seq;
}

sub Synonymous{
    my $line = $_[0];
    my ($CHROM, $POS, $REF, $NEW_POS, $ALT, $is_somatic, $MUT_TYPE, $STRAND, $NS_START, $PEP, $NUCLE, $PEP_ID, $ORF_start, $ORF_stop, $CONDON_START, $CONDON_END, $KMER_FREQ) = split /\t/, $line;
    my $front_end = $NEW_POS - $NS_START;
    my $back_start = length($ALT) + $front_end;
    my $front = substr $NUCLE, 0, $front_end;
    my $back = substr $NUCLE, $back_start;
    my $temp = $front.$REF.$back;
    my $orig_seq = NS2PEP($temp, $STRAND, length($PEP));
    if($orig_seq eq $PEP){ 
        return 1;
    }
    else{
        return 0;
    }
}

sub Writetmp{
    my ($headline, $out_f) = @_;
    open OUT, "> $out_f";
    print OUT $headline;
    foreach my $PEP(sort keys %remained){
        unless(exists $strings{$PEP}){
            print OUT $remained{$PEP}, "\n";
        }
    }
    close OUT;
}

sub RefProIndexing{
    my ($ref_f, $kmer) = @_;
    # my %Refs = ();
    my $seq;
    if($ref_f =~ /\.gz$/){ open REF, "gzip -dc $ref_f |" or die $!; }else{ open "< $ref_f" or die $!; }
    while(my $line=<REF>){
        chomp $line;
        if($line =~ /^>/){
            if($seq){ 
                my $kmer_num = length($seq) - $kmer + 1;
                for (my $n = 0; $n < $kmer_num; $n++) {
                    my $ind = substr $seq, $n, $kmer;
                    $Refs{$ind} = 1;
                }
                $seq = "";
            }
        }
        else{
            $seq .= $line;
        }
    }
    close REF;
}
