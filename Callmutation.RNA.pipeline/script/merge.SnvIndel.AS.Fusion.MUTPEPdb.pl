#!/usr/bin/perl -w
use strict;
die "\nUsage: $0   <SnvIndel_AS_Fusion.MUTPEP.list>   <output.unique.MUTPEP.fa>\n\n" unless (@ARGV == 2);

my %hash;
open (LIST, "$ARGV[0]") || die $!;
while (<LIST>) {
    chomp;
    open (IN, "$_") || die $!;
    $/ = ">"; <IN>;
    while (<IN>) {
        my @inf = split /\n/, $_;
        my $len = length($inf[1]);
        my $id = join("", (split /\_/, $inf[0])[1,2]);
        $hash{$len}{$inf[1]}{$id} = 1;
    }
    close IN;
    $/ = "\n";
}
close LIST;

open (OUT, "> $ARGV[1]") || die $!;
foreach my $len (sort{$a <=> $b} keys %hash) {
    foreach my $pep (sort{$a cmp $b} keys %{$hash{$len}}) {
        my $id = ">MUTPEP_" . join(".", sort{$a cmp $b} keys %{$hash{$len}{$pep}}) . "_len" . length($pep);
        print OUT "$id\n$pep\n";
    }
}
close OUT;
