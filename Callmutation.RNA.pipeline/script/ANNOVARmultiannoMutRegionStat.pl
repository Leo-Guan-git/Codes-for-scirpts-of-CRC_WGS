#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <multianno.txt>   <output>\n\n" unless (@ARGV == 2);

my %region;
my $all = 0;
open (IN, "$ARGV[0]") || die $!;
<IN>;
while (<IN>) {
    chomp;
    $all ++;
    my @inf = split /\t/, $_;
    $region{$inf[5]} ++;
}
close IN;

open (OUT, "> $ARGV[1]") || die $!;
print OUT join("\t", qw(Item  Number  Percentage.%))."\n";
foreach my $item (sort{$region{$b} <=> $region{$a}} keys %region) {
    my $rate = sprintf("%.2f", $region{$item}/$all*100);
    print OUT join("\t", $item, $region{$item}, $rate)."\n";
}
close OUT;
