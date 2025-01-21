#! /usr/bin/env perl -w
# @Author: Xiangyu Guan
# @Date:   2025-01-21
# @Last Modified by:   Xiangyu Guan
# @Last Modified time: 2025-01-21
use strict;
use warnings;
use POSIX;
use File::Basename;
use POSIX;
use List::Util;

@ARGV>=4 or die "Usage: perl $0 <Kmer file> <output file> <input file[s]>\n";

my ($output_f, @input_fs) = @ARGV;
my (%contents, %freq_sum);
my ($headline, $origin, $selected);
my $output_path = dirname $output_f;

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

sub Message{
    my $message = $_[0];
    my $datestring = localtime();
    print "\t $datestring >>> $message\n";
}

Message("Kmer Filter Job Start!");

system("mkdir -p $output_path") if (! -e $output_path);
if ($output_f =~ /\.gz$/){open OUT, "|gzip > $output_f" or die $!;}else{open OUT, ">$output_f" or die $!;}
my $input_f = $input_fs[0];
if ($input_f =~ /\.gz$/){open IN, "gzip -dc $input_f |" or die $!;}else{open IN, "<$input_f" or die $!;}
$headline = readline IN;
chomp $headline;
close IN;
print OUT "$headline\tPEP_MUT\tPEP_MUT_POS\n";

$origin = 0;
$selected = 0;
foreach my $input_f(@input_fs){
    Message("Input file: $input_f");
    if ($input_f =~ /\.gz$/){open IN, "gzip -dc $input_f |" or die $!;}else{open IN, "<$input_f" or die $!;}
    $headline = readline IN;
    chomp $headline;
    while (my $line=<IN>) {
        chomp $line;
        my ($CHROM, $POS, $REF, $NEW_POS, $ALT, $is_somatic, $MUT_TYPE, $STRAND, $NS_START, $PEP, $NUCLE, $PEP_ID, $ORF_start, $ORF_stop, $CONDON_START, $CONDON_END, $KMER_FREQ) = split /\t/, $line;
        $origin += 1;
        my ($Mut_aa_pos, $new_aa_seq, $raw_aa_seq, $PEP_MUT) = MutJuge($NEW_POS, $NS_START, $NUCLE, $REF, $ALT, $STRAND, $PEP);
        next if($PEP_MUT eq "Synonymous" || $Mut_aa_pos > length($PEP));
        # if($STRAND =~ /-/){
        #     $Mut_aa_pos = length($PEP) - ($Mut_aa_pos + length($new_aa_seq) -1) + 1;
        #     $new_aa_seq = scalar reverse $new_aa_seq;
        #     $raw_aa_seq = scalar reverse $raw_aa_seq;
        # }
        # if (exists $contents{$NUCLE}){ next;}
        # $contents{$NUCLE} = $_;
        $selected += 1;
        my @pep_mut_pos;
        if ($PEP_MUT =~ /^p\.DEL(\d+)([A-Z\*]+)$/){
            my $pre = $1 - 1;
            # my $after = $1 + length($raw_aa_seq);
            my $after = $1;
            $pep_mut_pos[0].= "$pre+$after";
        }
        else{
            for my $i(0..length($new_aa_seq)-1){
                push @pep_mut_pos, $Mut_aa_pos+$i;
            }
        }
        print OUT "$line\t$PEP_MUT\t".join(',', @pep_mut_pos)."\n";
        $freq_sum{$KMER_FREQ} += 1;
    }
    close IN;

}
Message("Got peptides number: $origin");
for my $freq(sort {$b <=> $a} keys %freq_sum){
    my $percentage = sprintf "%.4f", $freq_sum{$freq}/$origin*100;
    Message("Select $freq Kmer frequency peptides $freq_sum{$freq} ( $percentage\% )");
}

sub MutJuge{
    my ($NEW_POS, $NS_START, $NUCLE, $REF, $ALT, $STRAND, $PEP) = @_;
    my $Mut_ns_pos = $NEW_POS - $NS_START +1; #3
    my $Mut_aa_pos = ceil($Mut_ns_pos / 3); #1
    my ($raw_ns_seq, $new_ns_seq);
    if( (length($REF)-length($ALT)) % 3 != 0 ){
        $raw_ns_seq = substr($NUCLE, $Mut_aa_pos*3-3, $Mut_ns_pos-1-($Mut_aa_pos*3-3)).$REF.substr($NUCLE, $Mut_ns_pos+length($ALT)-1);
        $new_ns_seq = substr($NUCLE, $Mut_aa_pos*3-3, $Mut_ns_pos-1-($Mut_aa_pos*3-3)).$ALT.substr($NUCLE, $Mut_ns_pos+length($ALT)-1);
    }else{
        $raw_ns_seq = substr($NUCLE, $Mut_aa_pos*3-3, $Mut_ns_pos-1-($Mut_aa_pos*3-3)).$REF.substr($NUCLE, $Mut_ns_pos+length($ALT)-1, ceil(($Mut_ns_pos+length($REF)-1)/3)*3-1-($Mut_ns_pos+length($REF)-1)+1);
        $new_ns_seq = substr($NUCLE, $Mut_aa_pos*3-3, $Mut_ns_pos-1-($Mut_aa_pos*3-3)).$ALT.substr($NUCLE, $Mut_ns_pos+length($ALT)-1, ceil(($Mut_ns_pos+length($ALT)-1)/3)*3-1-($Mut_ns_pos+length($ALT)-1)+1);
    }
    # print "$raw_ns_seq\n$new_ns_seq\n";
    my ($new_aa_seq, $raw_aa_seq);
    my $num = -1;
    FIRST: for (my $i = 0; $i < List::Util::max((length($new_ns_seq), length($raw_ns_seq))); $i+=3) {
        my $n = $i / 3;
        my ($tmp_raw_ns_seq, $tmp_new_ns_seq) = (substr($raw_ns_seq, $i, 3), substr($new_ns_seq, $i, 3));
        if($STRAND =~ /-/){
            ($tmp_raw_ns_seq, $tmp_new_ns_seq) = (RverCom($tmp_raw_ns_seq), RverCom($tmp_new_ns_seq));
        }
        # my ($new_aa, $raw_aa) = ($AA_condon{substr($new_ns_seq, $i, 3)}, $AA_condon{substr($raw_ns_seq, $i, 3)});
        my ($new_aa, $raw_aa) = ($AA_condon{$tmp_new_ns_seq}, $AA_condon{$tmp_raw_ns_seq});
        if(($new_aa eq $raw_aa) && (($n - $num) == 1)){ 
            $Mut_aa_pos += 1;
            $num = $n;
            next FIRST;
        }else{
            # $new_aa_seq .= $AA_condon{substr($new_ns_seq, $i, 3)};
            # $raw_aa_seq .= $AA_condon{substr($raw_ns_seq, $i, 3)};
            $new_aa_seq .= $AA_condon{$tmp_new_ns_seq};
            $raw_aa_seq .= $AA_condon{$tmp_raw_ns_seq};
        }
    }
    if($STRAND =~ /-/){
        $Mut_aa_pos = length($PEP) - ($Mut_aa_pos + length($new_aa_seq) -1) + 1;
        $new_aa_seq = scalar reverse $new_aa_seq;
        $raw_aa_seq = scalar reverse $raw_aa_seq;
    }
    my $PEP_MUT;
    if($raw_aa_seq){
        # P.DEL4L                       DEL     缺失
        # P.4FETEFHSCC->LROSFTLVA       ->      移码
        # p.6A>V                        >       单氨基酸突变
        # P.IN512T                      INS     插入
        # p.7S*->FLK                    *       stop condon
        if($new_aa_seq){
            if(length($REF) == length($ALT) && length($REF) == 1){
                $PEP_MUT = "p.${Mut_aa_pos}${raw_aa_seq}>${new_aa_seq}";
            }
            else{
                $PEP_MUT = "p.${Mut_aa_pos}${raw_aa_seq}->${new_aa_seq}";
            }

        }
        else{
            $PEP_MUT = "p.DEL${Mut_aa_pos}${raw_aa_seq}";
        }
    }
    else{
        if($new_aa_seq){
            $PEP_MUT = "p.INS${Mut_aa_pos}${new_aa_seq}";
        }
        else{
            $PEP_MUT = "Synonymous"
        }
    }

    return($Mut_aa_pos, $new_aa_seq, $raw_aa_seq, $PEP_MUT);
}

sub RverCom{
    my $seq = $_[0];
    $seq =~ tr/ATCG/TAGC/;
    return(scalar reverse $seq);
}