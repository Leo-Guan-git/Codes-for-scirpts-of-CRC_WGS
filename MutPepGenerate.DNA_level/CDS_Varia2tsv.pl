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

my $scripts = basename $0;
@ARGV == 4 or die "Usage: perl $scripts <cds_from_genomic.fna[.gz]> <chromAlias.txt[.gz]> <mut_cds.faa[.gz]> <OutPut.tsv[.gz]>\n";

my ($CDS_fna_f, $chromAlias_f, $Mut_faa_f, $output_f) = @ARGV;
my (%chromAlias, %geneAlias);
if($chromAlias_f =~ /\.gz$/){ open ALIAS, "gzip -dc $chromAlias_f |" or die $!; }else{ open ALIAS, "< $chromAlias_f"; }
while (<ALIAS>) {
    if ($_ =~ /^#/){ next; }
    chomp;
    my @alias = split /\t/, $_;
    my ($UCSC, $Ensembl, $GenBank, $RefSeq) = split /\t/, $_;
    $chromAlias{$Ensembl} = $UCSC;
    $chromAlias{$GenBank} = $UCSC;
    $chromAlias{$RefSeq} = $UCSC;
    $chromAlias{$UCSC} = $UCSC;
}
close ALIAS;

if($CDS_fna_f =~ /\.gz$/){ open FNA, "gzip -dc $CDS_fna_f |" or die $!; }else{ open FNA, "< $CDS_fna_f" or die $!; }
while (<FNA>){
    next unless(/^>lcl\|([\w\.]+)\_cds\_([\w\.]+) \[gene=([^\]]+)\]/);
    my ($chrom_id, $cds_acccession, $gene_name) = ($1,$2,$3);
    $geneAlias{$cds_acccession} = $gene_name;
}
close FNA;

if($Mut_faa_f =~ /\.gz$/){ open IN, "gzip -dc $Mut_faa_f |" or die $!; }else{ open IN, "< $Mut_faa_f" or die $!; }
if($output_f =~ /\.gz$/){ open OUT, "| gzip -c > $output_f" or die $!; }else{ opem OUT, "> $output_f" or die $!; } 
print OUT "ID\tCHROM\tGene\tAA_POS\tAA_REF\tAA_ALT\tSTRAND\tNS_POS\tNS_REF\tNS_ALT\tInfo1\tInfo2\n";
while (<IN>) {
    chomp;
    next unless(/^>mut\|(([\w\.]+)\_prot\_([\w\.]+)\_mut\d+) \[strand: (\-?1)\] \[(\S+)\]$/);
    my ($mut_pro_id, $chrom, $cds_acccession, $strand) = ($1,$2,$3,$4);
    $chrom = $chromAlias{$chrom};
    my $gene = $geneAlias{$cds_acccession};
    if($strand == 1){
        $strand = "+";
    }else{
        $strand = "-";
    }
    for my $info(split /,/, $5){
        $info =~ /(\d+)([A-Z\*]+)>([A-Z\*]+)\(chr\w+\_(\d+)([ATCG]+)>([ATCG]+)\_([A-Z]+)\_([A-Z]+)\)+/;
        my ($AA_POS, $AA_REF, $AA_ALT, $NS_POS, $NS_REF, $NS_ALT, $mut_type1, $mut_type2) = ($1,$2,$3,$4,$5,$6,$7,$8);
        print OUT "$mut_pro_id\t$chrom\t$gene\t$AA_POS\t$AA_REF\t$AA_ALT\t$strand\t$NS_POS\t$NS_REF\t$NS_ALT\t$mut_type1\t$mut_type2\n";
    }
} 
close IN;
close OUT;
