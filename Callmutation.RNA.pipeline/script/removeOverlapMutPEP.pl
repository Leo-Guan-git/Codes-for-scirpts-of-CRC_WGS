#!/usr/bin/perl -w
use strict;
use File::Basename;

die "\nUsage: perl $0   <short.pep.list>   <long.pep.list>   <output>\n\n" unless (@ARGV == 3);

my @short;
open (AF, "$ARGV[0]") || die $!;
while (<AF>) {
    chomp;
    my @inf = split /\t/, $_;
    push @short, $inf[0];
}
close AF;

my @long;
open (BF, "$ARGV[1]") || die $!;
while (<BF>) {
    chomp;
    my @inf = split /\t/, $_;
    push @long, $inf[0];
}
close BF;

my $shortNum = scalar(@short);
my $longNum = scalar(@long);
my $shortFile = basename($ARGV[0]);
my $longFile = basename($ARGV[1]);
print "$shortFile vs $longFile  :  $shortNum vs $longNum\n";

my ($num1, $num2) = (0, 0);
open (OUT, "> $ARGV[2]") || die $!;
for (my $i=0; $i<@short; $i++) {
    my $flag = 0;
    for (my $j=0; $j<@long; $j++) {
        $flag ++ if ($long[$j] =~ /$short[$i]/);
    }
    print OUT "$short[$i]\n" if ($flag == 0);
    
    $num1 ++;
    $num2 ++;
    if ($num2 == 100) {
        my $date = `date`; chomp $date;
        print "processed $num1/$shortNum peptides  :  $date\n";
        $num2 = 0;
    }
}
close OUT;
