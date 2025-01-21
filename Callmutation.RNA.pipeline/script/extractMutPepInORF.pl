#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <mutant.peptides.gz>   <ORF.detail.info.gz>   <outdir>   <prefix>\n\n" unless (@ARGV == 4);

my ($outDir, $outPrefix) = @ARGV[2,3];

($ARGV[0] =~ /\.gz$/)?(open (MUT, "gzip -cd $ARGV[0] | ") || die $!):(open (MUT, "$ARGV[0]") || die $!);
my %head = ();
my @tmp = split /\t/, <MUT>; chomp $tmp[-1];
for (my $i=0; $i<@tmp; $i++) {
    $head{$tmp[$i]} = $i;
}
my @outHeader = @tmp;

my (%mutAA, %mutDetail);
while (<MUT>) {
    chomp;
    my @inf = split /\t/, $_;
    my $pep = $inf[$head{"mutAA"}];
    my $pepHead = substr($pep, 0, 4);
    $pepHead =~ s/I/L/g;
    $mutAA{$pepHead}{$pep} = 1;
    push @{$mutDetail{$pep}}, join("\t", @inf);
}
close MUT;


($ARGV[1] =~ /\.gz$/)?(open (ORF, "gzip -cd $ARGV[1] | ") || die $!):(open (ORF, "$ARGV[1]") || die $!);
%head = ();
@tmp = split /\t/, <ORF>; chomp $tmp[-1];
for (my $i=0; $i<@tmp; $i++) {
    $head{$tmp[$i]} = $i;
}
push @outHeader, @tmp;

open (OUT, " | gzip > $outDir/$outPrefix.SnvIndel.mutant.peptides.withORFmatch.gz") || die $!;
print OUT join("\t", @outHeader)."\n";

my ($block, $num) = (0, 0);
while (<ORF>) {
    chomp;
    my @inf = split /\t/, $_;
    my $aaSeq = $inf[$head{"aminoAcid"}];
    $aaSeq =~ s/I/L/g;
    foreach my $pepHead (sort{$a cmp $b} keys %mutAA) {
        next if ($aaSeq !~ /$pepHead/);
        foreach my $pep (sort{$a cmp $b} keys %{$mutAA{$pepHead}}) {
            my $transPEP = $pep;
            $transPEP =~ s/I/L/g;
            next if ($aaSeq !~ /$transPEP/);
            for (my $i=0; $i<@{$mutDetail{$pep}}; $i++) {
                print OUT join("\t", $mutDetail{$pep}[$i], @inf)."\n";
            }
        }
    }
    
    $num ++;
    $block ++;
    if ($block == 10000) {
        my $date = `date`; chomp $date;
        print "$inf[0]\t$num\t$date\n";
        $block = 0;
    }
}
close ORF;
close OUT;
