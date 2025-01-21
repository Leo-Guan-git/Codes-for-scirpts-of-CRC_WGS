#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <retained.peptides.list>   <mutant.peptides.noDBmatch.withORFmatch.unique.gz>   <output>\n\n" unless (@ARGV == 3);

my %target;
open (LIST, "$ARGV[0]") || die $!;
while (<LIST>) {
    chomp;
    my $file = $_;
    open (IN, "$file") || die $!;
    while (<IN>) {
        chomp;
        $target{$_} = 1;
    }
    close IN;
}
close LIST;

($ARGV[1] =~ /\.gz$/)?(open (BF, "gzip -cd $ARGV[1] | ") || die $!):(open (BF, "$ARGV[1]") || die $!);
($ARGV[2] =~ /\.gz$/)?(open (OUT, " | gzip > $ARGV[2]") || die $!):(open (OUT, " > $ARGV[2]") || die $!);
my @head = split /\t/, <BF>; chomp $head[-1];
print OUT join("\t", "mutPEPid", @head)."\n";
my $id = 0;
while (<BF>) {
    chomp;
    my @inf = split /\t/, $_;
    next if (!exists $target{$inf[0]});
    $id ++;
    print OUT join("\t", $id, @inf)."\n";
}
close BF;
close OUT;
