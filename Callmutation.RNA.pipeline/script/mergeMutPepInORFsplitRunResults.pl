#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <SnvIndel.mutant.peptides.withORFmatch.list>   <outdir>   <prefix>\n\n" unless (@ARGV == 3);

my ($outdir, $prefix) = @ARGV[1,2];

my (%mutSite, %gene, %region, %mutType, %mutDetail, %matchORF, %startCode),
open (LIST, "$ARGV[0]") || die $!;
open (OUT, " | gzip > $outdir/$prefix.withORFmatch.gz") || die $!;
my $fileNum = 0;
while (<LIST>) {
    chomp;
    $fileNum ++;
    open (IN, "gzip -cd $_ | ") || die $!;
    my @tmp = split /\t/, <IN>; chomp $tmp[-1];
    print OUT join("\t", @tmp)."\n" if ($fileNum == 1);
    my %head;
    for (my $i=0; $i<@tmp; $i++) {
        $head{$tmp[$i]} = $i;
    }
    while (<IN>) {
        chomp;
        my @inf = split /\t/, $_;
        print OUT join("\t", @inf)."\n";
        my $mutAA = $inf[$head{"mutAA"}];
        my $site = join(".", @inf[$head{"Chr"}, $head{"Start"}, $head{"End"}]);
        $mutSite{$mutAA}{$site} = 1;
        
        my $ORF = join(".", @inf[$head{"transID"}, $head{"transOrfID"}, $head{"strand"}, $head{"frame"}, $head{"start"},$head{"end"}]);
        $matchORF{$mutAA}{$ORF} = 1;
        
        $gene{$mutAA}{$inf[$head{"Gene"}]} = 1;
        $region{$mutAA}{$inf[$head{"Region"}]} = 1;
        $mutType{$mutAA}{$inf[$head{"mutType"}]} = 1;
        $mutDetail{$mutAA}{$inf[$head{"mutDetail"}]} = 1;
        $startCode{$mutAA}{$inf[$head{"startCode"}]} = 1;
    }
    close IN;
}
close LIST;

open (OUT, " | gzip > $outdir/$prefix.withORFmatch.unique.gz") || die $!;
print OUT join("\t", qw(mutAA  mutAAlength  mutSiteNum  mutSites  geneNum  genes  regionNum  regions  mutTypeNum  mutTypes  mutDetailNum  mutDetails  startCodeNum  startCodes  matchORFnum  matchORFs))."\n";

foreach my $mutAA (sort{$a cmp $b} keys %mutSite) {
    my $num = scalar(keys %{$mutSite{$mutAA}});
    my $detail = join("|", sort{$a cmp $b} keys %{$mutSite{$mutAA}});
    my @outLine = ($num, $detail);
    $num = scalar(keys %{$gene{$mutAA}});
    $detail = join("|", sort{$a cmp $b} keys %{$gene{$mutAA}});
    push @outLine, $num, $detail;
    
    $num = scalar(keys %{$region{$mutAA}});
    $detail = join("|", sort{$a cmp $b} keys %{$region{$mutAA}});
    push @outLine, $num, $detail;
    
    $num = scalar(keys %{$mutType{$mutAA}});
    $detail = join("|", sort{$a cmp $b} keys %{$mutType{$mutAA}});
    push @outLine, $num, $detail;
    
    $num = scalar(keys %{$mutDetail{$mutAA}});
    $detail = join("|", sort{$a cmp $b} keys %{$mutDetail{$mutAA}});
    push @outLine, $num, $detail;
    
    $num = scalar(keys %{$startCode{$mutAA}});
    $detail = join("|", sort{$a cmp $b} keys %{$startCode{$mutAA}});
    push @outLine, $num, $detail;
    
    $num = scalar(keys %{$matchORF{$mutAA}});
    $detail = join("|", sort{$a cmp $b} keys %{$matchORF{$mutAA}});
    push @outLine, $num, $detail;
    
    print OUT join("\t", $mutAA, length($mutAA), @outLine)."\n";
}
close OUT;