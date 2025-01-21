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
use Fcntl qw( SEEK_SET );
use Getopt::Std;
# use Data::Dumper;
use vars qw($opt_s $opt_i $opt_g $opt_b $opt_c $opt_l $opt_p $opt_a);
getopts('s:i:g:b:c:l:p:a:');
my $version = "[v1.0]\nGenerate mutate protein from vcf and cds files";
# my $fa_fai_f = $fa_f.".fai";
my ($snv_f, $indel_f, $snp_f, $bin_size, $cds_faa_f, $out_prefix, $mut_out_f, $cds_mut_out_f, $pep_mut_out_f, $LINEBASES, $log_out_f, $message, $chromAlias_f);
my (%chromAlias, %variations);
$bin_size = 100000;
$LINEBASES = 80;
$out_prefix = "./test.mut";
$chromAlias_f = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/hg38/gtf/hg38.p13.chromAlias.txt";
$cds_faa_f = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/hg38/fasta/pep/GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna.gz";
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

if(! ($opt_s || $opt_i || $opt_g)){
    print "Usage: $0 $version\n";
    print "--------------------------------------------------------------------\n";
    print "-s somatic snv file..................................<snv.vcf[.gz]>\n";
    print "-i somatic indel file................................<indel.vcf[.gz]>\n";
    print "-g snp file..........................................<*.vcf[.gz]>\n";
    print "-b bin size for chromosme region segementation........[100000]\n";
    print "-c coding region sequence file........................[GRCh38.p13_cds_from_genomic.fna.gz]\n";
    print "-l line bases for fasta file output...................[80]\n";
    print "-p Prefix.............................................[./test.mut]\n";
    print "-a chromosme id alias file............................[hg38.p13.chromAlias.txt]\n";
    die "Please Check your parameters\n";    
}

# $snv_f = $opt_s if($opt_s);
# $indel_f = $opt_i if($opt_i);
# $snp_f = $opt_g if($opt_g);
$bin_size = $opt_b if($opt_b);
$cds_faa_f = $opt_c if($opt_c);
$out_prefix = $opt_p if($opt_p);
$LINEBASES = $opt_l if($opt_l);
$chromAlias_f = $opt_a if($opt_a);
$mut_out_f = $out_prefix.".all.vcf.gz";
$cds_mut_out_f = $out_prefix.".cds.fna.gz";
$pep_mut_out_f = $out_prefix.".cds.faa.gz";
$log_out_f = $out_prefix.".log";
open LOG, "> $log_out_f";

if($opt_s){
    $snv_f = $opt_s;
    if (! -e $snv_f){
        Message("somatic snv file $snv_f doesn't not exist.");
        die "somatic snv file $snv_f doesn't not exist.";
    }
}

if($opt_i){
    $indel_f = $opt_i;
    if (! -e $indel_f){
        $message = Message("somatic indel file $indel_f doesn't not exist.");
        die "somatic indel file $indel_f doesn't not exist.";
    }
}

if($opt_g){
    $snp_f = $opt_g;
    if (! -e $snp_f){
        $message = Message("germline snp file $snp_f doesn't not exist.");
        die "germline snp file $snp_f doesn't not exist.";
    }
}

=ignored
# generate fai file if not exist
unless(-e $fa_fai_f){
    # Message("Generate $fa_fai_f file");
    my $samtools = `which samtools`;
    $samtools = "/share/app/samtools/1.11/bin/samtools" unless($samtools);
    # Message("Command Line: $samtools faidx $fa_f");
    system("$samtools faidx $fa_f");
    # Message("samtools index generated")
}
=cut
Message("chromosme Alias file $chromAlias_f parsing start...");
# my %chromAlias;
open ALIAS, "< $chromAlias_f";
while (<ALIAS>) {
    if ($_ =~ /^#/){ next; }
    chomp;
    my @alias = split /\t/, $_;
    my ($UCSC, $Ensembl, $GenBank, $RefSeq) = split /\t/, $_;
    $chromAlias{$Ensembl} = $UCSC;
    $chromAlias{$GenBank} = $UCSC;
    $chromAlias{$RefSeq} = $UCSC;
    $chromAlias{$UCSC} = $UCSC
}
close ALIAS;
Message("chromosme Alias file $chromAlias_f parsing done");
=ignored
my %faidx;
open FAI, "<$fa_fai_f";
while (<FAI>){
    chomp $_;
    my ($CHROM, $LENGTH, $OFFSET, $LINEBASES, $LINEWIDTH) = split /\t/, $_;
    $faidx{$CHROM}{"LENGTH"} = $LENGTH; 
    $faidx{$CHROM}{"OFFSET"} = $OFFSET; 
    $faidx{$CHROM}{"LINEBASES"} = $LINEBASES; 
    $faidx{$CHROM}{"LINEWIDTH"} = $LINEWIDTH;
}
close FAI;
=cut

if($opt_s){
    Message("input file $snv_f parsing start...");
    my $num = VcfFileRead($snv_f, $bin_size, "SOMATIC");
    Message("input file $snv_f parsing done.");
    Message("$num positions are selected");
}
if($opt_i){
    Message("input file $indel_f parsing start...");
    my $num = VcfFileRead($indel_f, $bin_size, "SOMATIC");
    Message("input file $indel_f parsing done.");
    Message("$num positions are selected");
}
if($opt_g){
    Message("input file $snp_f parsing start...");
    my $num = VcfFileRead($snp_f, $bin_size, "GERMLINE");
    Message("input file $snp_f parsing done.");
    Message("$num positions are selected");
}

Message("All vcf info write to $mut_out_f...");
if ($mut_out_f =~ /\.gz$/){ open MUT, "| gzip -c > $mut_out_f" or die $!; }else{ open MUT, ">$mut_out_f" or die $!; }
print MUT "#CHROM\tPOS\tREF\tALT\tTYPE1\tTYPE2\n";
foreach my $chrom(sort keys %variations){
    foreach my $bin(sort {$a <=> $b} keys %{$variations{$chrom}}){
        foreach my $pos(sort {$a <=> $b} keys %{$variations{$chrom}{$bin}}){
            print MUT "$variations{$chrom}{$bin}{$pos}\n";
        }
    }
}
close MUT;
Message("vcf info write done");

Message("start CDS faa file Parsing...");
if($cds_faa_f =~ /\.gz$/){ open FAA, "gzip -dc $cds_faa_f | " or die $!; }else{ open FAA, "< $cds_faa_f" or die $!; }
if($cds_mut_out_f =~ /\.gz$/){ open CDS, " | gzip -c > $cds_mut_out_f" or die $!; }else{ open CDS, "> $cds_mut_out_f" or die $!;}
if($pep_mut_out_f =~ /\.gz$/){ open PEP, " | gzip -c > $pep_mut_out_f" or die $!; }else{ open PEP, "> $pep_mut_out_f" or die $!;}
my $mut_num = 0;
my $pep_num = 0;
my ($chrom, $cds, $seq, $strand);
my @LOCATIONS;
while (my $line=<FAA>) {
    chomp $line;
    if ($line =~ /^>lcl\|([\w\.]+)\_cds\_([\w\.]+) /){
        # print $_, "\n";
        $pep_num += 1;
        Message("$pep_num proteins are parsed") if($pep_num % 10000 ==0);
        if ($seq){
            # print length($seq);
            # print "$chrom\t$cds\t$mut_num\n$seq\n@LOCATIONS\n";
            my ($new_seq, @NS_Mutations) = MutSeqGenerate($chrom, $seq, $strand, @LOCATIONS);
            if ($#NS_Mutations > -1){
                my ($mut_pep, @PEP_Mutations) = MutPepGenerate($strand, $new_seq, @NS_Mutations);
                if ($mut_pep && ($#PEP_Mutations > -1)){
                    $mut_num +=1;
                    print CDS ">mut|$chrom\_prot\_$cds\_mut$mut_num [".join(",", @NS_Mutations)."]\n";
                    print CDS SeqLineSplit($new_seq, $LINEBASES), "\n";
                    print PEP ">mut|$chrom\_prot\_$cds\_mut$mut_num [strand: $strand] [".join(",", @PEP_Mutations)."]\n";
                    print PEP SeqLineSplit($mut_pep, $LINEBASES), "\n";
                    Message("$mut_num proteins are generated") if($mut_num % 10000 ==0);
                }
            }
        }
        ($chrom, $cds) = ($1,$2);
        $seq = "";
        @LOCATIONS = ();
        $strand = 1;
        $line =~ /\[location=(complement\()?(join\()?([\d\.\,\>\<]+)\)?\)?\]/;
        my $location_info = $3;
        if ($1){ $strand = -1;}
        $location_info =~ tr/\<\>//d;
        my @locations = split /,/, $location_info;
        foreach my $loc(@locations){
            my ($start, $end);
            if ($loc =~ /\.{2}/){ ($start, $end) = split /\.\./, $loc; } else { ($start, $end) = ($loc, $loc); }
            push @LOCATIONS, [$start, $end];
        }
    }else{
        $seq .= $line;
    }
}
if ($seq){
    my ($new_seq, @NS_Mutations) = MutSeqGenerate($chrom, $seq, $strand, @LOCATIONS);
    if ($#NS_Mutations > -1){
        my ($mut_pep, @PEP_Mutations) = MutPepGenerate($strand, $new_seq, @NS_Mutations);
        if ($mut_pep && ($#PEP_Mutations > -1)){
            $mut_num +=1;
            print CDS ">mut|$chrom\_prot\_$cds\_mut$mut_num [".join(",", @NS_Mutations)."]\n";
            print CDS SeqLineSplit($new_seq, $LINEBASES), "\n";
            print PEP ">mut|$chrom\_prot\_$cds\_mut$mut_num [strand: $strand] [".join(",", @PEP_Mutations)."]\n";
            print PEP SeqLineSplit($mut_pep, $LINEBASES), "\n";
            Message("$mut_num proteins are generated") if($mut_num % 10000 ==0);
        }
    }
}
close FAA;
close CDS;
close PEP;
Message("CDS faa file Parse and Mut protein generate are done");
Message("$mut_num mutated proteins are generated from $pep_num origin proteins");
=ingored
sub FastaSeqRead{
    # in order to fit for complete human genome.
    # do not load a big reference to memory.
    my ($fa_f,$CHROM) = @_;
    if($fa_f =~ /\.gz$/){ open REF, "gzip -dc $fa_f |" or die $!; }else{ open REF, "< $fa_f" or die $!; }
    my $seq;
    my $start_byte = $faidx{$CHROM}{"OFFSET"};
    my $length_byte = floor($faidx{$CHROM}{"LENGTH"} / $faidx{$CHROM}{"LINEBASES"}) + $faidx{$CHROM}{"LENGTH"};
    seek REF, $start_byte, SEEK_SET or die $!;
    read REF, $seq, $length_byte;
    $seq =~ tr/\n//d;
    $seq =~ tr/a-z/A-Z/;
    close REF;
    return $seq
}
=cut
close LOG;
sub Message{
    my $message = $_[0];
    my $datestring = localtime();
    $message = "$datestring >>> $message\n";
    print $message;
    print LOG $message;
}

sub VcfFileRead{
    # select [snvs] and [non-frameshift indels] to generate personalized
    # CDS and Protein sequences.

    my ($vcf_f, $bin_size, $mut_type) = @_;
    my $desc;
    my $num = 0;
    if($vcf_f =~ /\.gz$/){ open VCF, "gzip -dc $vcf_f |" or die $!; }else{ open VCF, "<$vcf_f" or die $!; }
    while (<VCF>){
        if (/^#/){next;}
        chomp;
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO) = (split /\t/, $_)[0..7];
        if ($FILTER ne "PASS"){ next; }
        my $bin_num = floor $POS/$bin_size;
        if (exists $variations{$CHROM}{$bin_num}{$POS}){ next; }
        if ($ALT =~ /,/){ $ALT = (split /,/, $ALT)[0]; }
        if ((length($ALT)-length($REF)) % 3 != 0) { next; }
        $variations{$CHROM}{$bin_num}{$POS} = "$CHROM\t$POS\t$REF\t$ALT";
        if ($mut_type eq "GERMLINE"){
            $variations{$CHROM}{$bin_num}{$POS} .= "\tGERMLINE";
        }else{
            $variations{$CHROM}{$bin_num}{$POS} .= "\tSOMATIC";
        }
        if (length($ALT)==length($REF)){
            $variations{$CHROM}{$bin_num}{$POS} .= "\tSNV";
        }elsif(length($ALT)>length($REF)){
            $variations{$CHROM}{$bin_num}{$POS} .= "\tINS";
        }else{
            $variations{$CHROM}{$bin_num}{$POS} .= "\tDEL";
        }
        $num += 1;
    }
    close VCF;
    return $num;
}

sub NS2PEP{
    my ($ns_seq, $strand) = @_;
    if (length($ns_seq) % 3 != 0){return "";}
    my $pep_seq;
    # $ns_seq =~ tr/atcg/ATCG/;
    if ($strand == -1){
        $ns_seq =~ tr/ATCG/TAGC/;
        $ns_seq = reverse $ns_seq;
    }
    for (my $n = 0; $n < length($ns_seq); $n += 3) {
        my $condon = substr $ns_seq, $n, 3;
        # if (! exists $AA_condon{$condon}){print $condon,"\n";}
        $pep_seq .= $AA_condon{$condon};
    }
    return $pep_seq;
}

sub MutSeqGenerate{
    my ($chrom, $seq, $strand, @LOCATIONS) = @_;
    $seq =~ tr/atcg/ATCG/;
    if ($strand == -1){
        $seq =~ tr/ATCG/TAGC/;
        $seq = reverse $seq;
    }
    my (@tmp_locations, @new_locations, @NS_Mutations);
    my ($tmp_seq, $tmp_length, $new_seq, $desc);
    my $pointer = 0;
    my $new_pos = 0;
    # push @new_locations, $new_pos;
    for (my $n = 0; $n < $#LOCATIONS+1; $n++) {
        my ($start, $end) = ($LOCATIONS[$n][0], $LOCATIONS[$n][1]);
        $tmp_seq = "";
        my $pointer_rela = 0;
        # @tmp_locations = ();
        # push @tmp_locations, $pointer_rela;
        my $bin_min = floor $start/$bin_size;
        my $bin_max = floor $end/$bin_size;
        # push @tmp_locations, $start;
        OUTER: foreach my $bin($bin_min..$bin_max) {
            foreach my $pos(sort {$a <=> $b} keys %{$variations{$chromAlias{$chrom}}{$bin}}){
                if ($pos < $start) { next; }
                if ($pos > $end) { last OUTER; }
                $tmp_length = $pos-$start - $pointer_rela;
                # $tmp_length = $pos - $tmp_locations[$#tmp_locations] + 1;
                $tmp_seq .= substr $seq, $pointer, $tmp_length;
                $pointer_rela += $tmp_length;
                $pointer += $tmp_length;
                $new_pos += $tmp_length;
                # push @tmp_locations, $pointer;
                # push @new_locations, $new_pos;
                my ($chrom, $pos, $ref, $alt, $type1, $type2) = split /\t/, $variations{$chromAlias{$chrom}}{$bin}{$pos};
                $tmp_seq .= $alt;
                my $desc = $new_pos+1;
                $desc .= "\_$ref>$alt($chrom\_$pos\_$type1\_$type2)";
                $pointer += length($ref);
                $pointer_rela += length($ref);
                $new_pos += length($alt);
                # push @tmp_locations, $pointer;
                # push @new_locations, $new_pos;
                push @NS_Mutations, $desc;
            }
        }
        if (($end-$start+1) > $pointer_rela){
            $tmp_length = $end - $start - $pointer_rela + 1;
            $tmp_seq .= substr $seq, $pointer, $tmp_length;
            $pointer += $tmp_length;
            $new_pos += $tmp_length;
            # push @tmp_locations, $pointer;
            # push @new_locations, $new_pos;
        }
        $new_seq .= $tmp_seq;
    }
    return($new_seq, @NS_Mutations);
}

sub Locate{
    my ($from, $to) = @_;
    my @juge;
    push @juge, [$from, $from-1, $from-2];
    push @juge, [$to+2, $to+1, $to];
    return($juge[0][$from % 3], $juge[1][$to % 3]) 
}

sub MutPepGenerate{
    my ($strand, $new_seq, @NS_Mutations) = @_;
    my $new_pep = NS2PEP($new_seq, $strand);
    if (! $new_pep){ return("",""); }
    if ($new_pep =~ /\*$/){$new_pep = substr $new_pep, 0, -1;}
    if ($new_pep =~ /\*/){ return("",""); }
    my @PEP_Mutations;
    foreach (@NS_Mutations){
        $_ =~ /(\d+)_([ATCG]+)>([ATCG]+)\((\w+)_(\d+)_(SOMATIC|GERMLINE)_(SNV|INS|DEL)\)/;
        my ($chr, $ns_pos, $cds_pos, $ref, $alt, $type1, $type2) = ($4, $5, $1, $2, $3, $6, $7);
        my $ns_from = $cds_pos-1;
        my $ref_to = $ns_from+length($ref)-1;
        my $alt_to = $ns_from+length($alt)-1;
        my ($ref_ind1, $ref_ind2) = Locate($ns_from, $ref_to);
        # print "$ref_ind1\t$ref_ind2\n";
        my ($alt_ind1, $alt_ind2) = Locate($ns_from, $alt_to);
        # print "$alt_ind1\t$alt_ind2\n";
        my $tmp_seq = substr($new_seq, $ref_ind1, $ns_from-$ref_ind1).$ref.substr($new_seq, $ns_from+length($alt), $ref_ind2-$ref_to);
        my $ref_pep = NS2PEP($tmp_seq, $strand);
        my $alt_pep = NS2PEP(substr($new_seq, $alt_ind1, $alt_ind2-$alt_ind1+1), $strand);
        if ($ref_pep eq $alt_pep){ next; }
        my $ind;
        if ($strand == 1){$ind = floor($alt_ind1/3) + 1;}
        if ($strand == -1){$ind = length($new_pep) - floor($alt_ind2/3) + 1;} 
        push @PEP_Mutations, "$ind$ref_pep>$alt_pep($chr\_$ns_pos$ref>$alt\_$type1\_$type2)";
    }
    return($new_pep, @PEP_Mutations);
}

sub SeqLineSplit{
    my ($seq, $LINEBASES) = @_;
    my @seq_splited;
    for (my $n = 0; $n < length($seq); $n += $LINEBASES) {
        push @seq_splited, substr($seq, $n, $LINEBASES);

    }
    return join("\n", @seq_splited);
}