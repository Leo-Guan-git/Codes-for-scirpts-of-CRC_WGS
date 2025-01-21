#! /usr/bin/env perl -w
# @Author: Xiangyu Guan
# @Date:   2025-01-21
# @Last Modified by:   Xiangyu Guan
# @Last Modified time: 2025-01-21

=pod
This pepline is used to process mutate peptides that are filtered
through ORF prediction.
output peptides meet 2 conditions:
    1. split to multi-files by peptide length;
    2. deduplicte by their nucleotide sequence in each length;
    3. filter pseudo sequences that do not supported in Original data
=cut

use strict;
use warnings;
use POSIX;
use File::Basename;
use Data::Dumper;

@ARGV==3 or die "Usage: perl $0 <output file Prefix> <input file>\n";

my ($output_prefix, $bam_f, $input) = @ARGV;
my (%tmp, %contents, %strings, %remained, %origin, %output_fh, %pep_bin, %Kmer_bin, %selected_reads);
my ($headline, $pep_length, $output_path);
my ($samtools, $chrom_pos);
# my @selected_reads;

if(-e "/share/app/samtools/1.11/bin/samtools"){
    $samtools = "/share/app/samtools/1.11/bin/samtools";
}

Message("step2 filter Job Start!");
Message("Input file: $input");
$output_path = dirname $output_prefix;
if ( ! -e $output_path){system "mkdir -p $output_path";}
if ($input =~ /\.gz$/){open IN, "gzip -dc $input |" or die $!;}else{open IN, "< $input";}
$headline = readline IN;
chomp $headline;
$chrom_pos = 0;
while (my $line=<IN>){
    chomp $line;
    my ($CHROM, $POS, $REF, $NEW_POS, $ALT, $is_somatic, $MUT_TYPE, $STRAND, $NS_START, $PEP, $NUCLE, $PEP_ID, $ORF_start, $ORF_stop, $CONDON_START, $CONDON_END) = split /\t/, $line;
    $pep_length = length $PEP;
    if (! exists $tmp{$pep_length}){
        open $output_fh{$pep_length}, "> $output_prefix.PEP$pep_length.somatic.tsv" or die $!;
        print { $output_fh{$pep_length} } "$headline\tKMER_FREQ\n";
        $tmp{$pep_length} = 1;
    }
    # if (! exists $tmp{$pep_length}){ $tmp{$pep_length} = 1;}
    # if (exists $contents{$NUCLE}){next;}
    if ($is_somatic==0){next;}
    if ($chrom_pos eq "$CHROM\t$POS") {
        $pep_bin{$pep_length}{$NUCLE} = $line;
    }else{
        if ($chrom_pos) {
            my ($chrom, $start) = split /\t/, $chrom_pos;
            # @selected_reads = ExtractReads($bam_f, $chrom, $start);
            ExtractReads($bam_f, $chrom, $start);
            # %Kmer_bin = ();
            # undef %Kmer_bin;
            GenerateKmer(%pep_bin);
            for my $len(sort {$a <=> $b} keys %pep_bin){
                for my $kmer(sort {$a cmp $b} keys %{$pep_bin{$len}}){
                    if(exists $Kmer_bin{$kmer}){
                        # for my $record(@{$pep_bin{$len}{$kmer}}){
                            # print { $output_fh{$len} } "$record\t$Kmer_bin{$kmer}\n";
                        # }
                        print { $output_fh{$len} } "$pep_bin{$len}{$kmer}\t$Kmer_bin{$kmer}\n";
                        $remained{$len} += 1;
                    }
                }
            }
        }
        $chrom_pos = "$CHROM\t$POS";
        # %pep_bin = ();
        undef %pep_bin;
        $pep_bin{$pep_length}{$NUCLE} = $line;
    }
    # print { $output_fh{$pep_length} } "$line\n";
    # $contents{$NUCLE} = 1;
    # $remained{$pep_length} += 1;
}

if ($chrom_pos) {
    my ($chrom, $start) = split /\t/, $chrom_pos;
    # @selected_reads = ExtractReads($bam_f, $chrom, $start);
    ExtractReads($bam_f, $chrom, $start);
    %Kmer_bin = ();
    GenerateKmer(%pep_bin);
    for my $len(sort {$a <=> $b} keys %pep_bin){
        for my $kmer(sort {$a cmp $b} keys %{$pep_bin{$len}}){
            if(exists $Kmer_bin{$kmer}){
                # for my $record(@{$pep_bin{$len}{$kmer}}){
                    # print { $output_fh{$pep_length} } "$record\t$Kmer_bin{$kmer}\n";
                # }
                print { $output_fh{$len} } "$pep_bin{$len}{$kmer}\t$Kmer_bin{$kmer}\n";
                $remained{$len} += 1;
            }
        }
    }
}

foreach my $out_fh(keys %output_fh){close $output_fh{$out_fh};}
foreach my $len(sort {$a <=> $b} keys %remained){ Message("got $remained{$len} $len length peptide from input file"); }

sub Message{
    my $message = $_[0];
    my $datestring = localtime();
    print "\t $datestring >>> $message\n";
}

sub ExtractReads{
    my ($bam_f, $chr, $start) = @_;
    undef %selected_reads;
    # my %Reads;
    open BAM, "$samtools view $bam_f $chr:$start-$start |" or die $!;
    while (<BAM>) {
        chomp;
        my $read = (split /\t/,)[9];
        $selected_reads{$read} = 1;
    }
    close BAM;
    # my @selected_reads = keys %Reads;
    # undef %Reads;
    # return(@selected_reads);
}

sub GenerateKmer{
    my (%hash) = @_;
    undef %Kmer_bin;
    for my $pep_length(sort {$a <=> $b} keys %hash){
        my $nucle_length = 3 * $pep_length;
        # for my $read(@selected_reads){
        for my $read(keys %selected_reads){
            my $num = length($read) - $nucle_length + 1;
            for (my $i = 0; $i <= $num; $i++) {
                my $kmer = substr $read, $i, $nucle_length;
                $Kmer_bin{$kmer} += 1;
            }
        }
    }
}
