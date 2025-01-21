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

@ARGV >= 5 or die "Usage: perl $0 <annovar.vcf> <sample id> <pep.tsv> <Outdir> <protocol> [<genename.list>]\n";
my ($ann_vcf_f, $sampid, $pep_tsv_f, $Outdir, $protocol, $gene_list_f) = @ARGV;
my (%Genes, %Annos);

if ($gene_list_f){
    open IN, "< $gene_list_f" or die $!;
    while (<IN>) {
        chomp;
        my ($CHROM, $GENE_ID, $GENE_NAME) = split /\t/;
        $Genes{$CHROM}{$GENE_ID} = $GENE_NAME;
    }
    close IN
}

if ($ann_vcf_f =~ /\.gz$/){ open ANN, "gzip -dc $ann_vcf_f |" or die $!; }else{ open ANN, "< $ann_vcf_f" or die $!; }
while (<ANN>) {
    chomp;
    next if(/^#/);
    my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO) = split /\t/;
    $INFO =~ /ANNOVAR_DATE=[^;]+;Func\.\Q$protocol\E=([^;]+);Gene\.\Q$protocol\E=([^;]+);GeneDetail\.\Q$protocol\E=([^;]+);ExonicFunc\.\Q$protocol\E=([^;]+);AAChange\.\Q$protocol\E=([^;]+);ALLELE_END$/;
    my ($Func, $Gene, $GeneDetail, $ExonicFunc, $AAChange) = ($1,$2,$3,$4,$5);
    $Func =~ s/\\x3b/;/g;
    $Gene =~ s/\\x3b/;/g;
    $GeneDetail =~ s/\\x3b/;/g;
    $GeneDetail =~ s/\\x3d/=/g;
    $ExonicFunc =~ s/\\x3b/;/g;
    $AAChange =~ s/\\x3b/;/g;
    my $Gene_name = "-";
    if ($gene_list_f){
        my @Gene_name;
        for my $gene(split /;/, $Gene){
            push @Gene_name, $Genes{$CHROM}{$gene};
        }
        $Gene_name = join ";", @Gene_name;
    }
    $Annos{"$CHROM\_$POS\_$REF\_$ALT"}{Func} = $Func;
    $Annos{"$CHROM\_$POS\_$REF\_$ALT"}{Gene_id} = $Gene;
    $Annos{"$CHROM\_$POS\_$REF\_$ALT"}{Gene_name} = $Gene_name;
    $Annos{"$CHROM\_$POS\_$REF\_$ALT"}{GeneDetail} = $GeneDetail;
    $Annos{"$CHROM\_$POS\_$REF\_$ALT"}{ExonicFunc} = $ExonicFunc;
    $Annos{"$CHROM\_$POS\_$REF\_$ALT"}{AAChange} = $AAChange;
}
close ANN;
my %ANNOVA_precedence = ();
if ($pep_tsv_f =~ /\.gz$/){ open TSV, "gzip -dc $pep_tsv_f |" or die $!; }else{ open TSV, "< $pep_tsv_f" or die $!; }
my $outfile = "$Outdir/$sampid.PEP.final.annot.tsv";
open OUT, "> $outfile" or die $!;
my $headline= readline TSV;
$headline =~ s/[\n\r]+$//;
print OUT "$headline\tFunc\tGene_id\tGene_name\tGeneDetail\tExonicFunc\tAAChange\n";
while (my $line=<TSV>) {
    chomp $line;
    # my ($CHROM, $POS, $REF, $ALT) = (split /\t/, $line)[(2,3,4,6)];
    my ($CHROM, $POS, $REF, $ALT) = (split /\t/, $line)[(1,2,3,5)];
    my $ind = "$CHROM\_$POS\_$REF\_$ALT";
    print OUT "$line\t$Annos{$ind}{Func}\t$Annos{$ind}{Gene_id}\t$Annos{$ind}{Gene_name}\t$Annos{$ind}{GeneDetail}\t$Annos{$ind}{ExonicFunc}\t$Annos{$ind}{AAChange}\n";
}
close TSV;
close OUT;
