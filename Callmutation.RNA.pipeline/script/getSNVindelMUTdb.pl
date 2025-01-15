#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <SnvIndel.mutant.peptides.noDBmatch.withORFmatch.unique.noOverlap.gz>   <outdir>   <outprefix>\n\n" unless (@ARGV == 3);

my ($outdir, $prefix) = @ARGV[1,2];

open (IN, "gzip -cd $ARGV[0] | ") || die $!;
my @tmp = split /\t/, <IN>; chomp $tmp[-1];
my %head;
for (my $i=0; $i<@tmp; $i++) {
    $head{$tmp[$i]} = $i;
}

my (%pepMutType, %pepFeature);
while (<IN>) {
    chomp;
    my @inf = split /\t/, $_;
    my ($pep, $len) = @inf[$head{"mutAA"}, $head{"mutAAlength"}];
    
    my @types = split /\|/, $inf[$head{"mutTypes"}];
    for (my $i=0; $i<@types; $i++) {
        my $tag = (split /\./, $types[$i])[-1];
        $tag = uc($tag);
        push @{$pepMutType{$len}{$pep}}, $tag;
    }
    my ($mutSites, $genes, $regions, $mutTypes, $mutDetails, $startCodes, $matchORFs) = @inf[$head{"mutSites"}, $head{"genes"}, $head{"regions"}, $head{"mutTypes"}, $head{"mutDetails"}, $head{"startCodes"}, $head{"matchORFs"}];
    $pepFeature{$pep} = join("\t", $mutSites, $genes, $regions, $mutTypes, $mutDetails, $startCodes, $matchORFs);
}
close IN;

open (OUT, "> $outdir/$prefix.SnvIndel.MUTPEP.fasta") || die $!;
open (OUT2, "> $outdir/$prefix.SnvIndel.MUTPEP.features.tsv") || die $!;
print OUT2 join("\t", qw(MUTPEPid  MUTPEPseq  mutSites  genes  regions  mutTypes  mutDetails  startCodes  matchORFs))."\n";
my $idNum = 0;
foreach my $len (sort{$a <=> $b} keys %pepMutType) {
    foreach my $pep (sort{$a cmp $b} keys %{$pepMutType{$len}}) {
        $idNum ++;
        my $mutTag = join(".", sort{$a cmp $b} @{$pepMutType{$len}{$pep}});
        my $faID = join("_", ">MUTPEP", $mutTag, $idNum);
        print OUT "$faID\n$pep\n";
        
        $faID =~ s/^>//;
        print OUT2 join("\t", $faID, $pep, $pepFeature{$pep})."\n";
    }
}
close OUT;
close OUT2;
