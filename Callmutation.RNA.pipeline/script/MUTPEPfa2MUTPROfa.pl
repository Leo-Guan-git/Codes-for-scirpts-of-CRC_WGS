#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <INPUT(MUTPEP.fa)>   <outMUTPROfa>   <outMUTPRO2MUTPEPpos>   [refPROdb.fa]\n\n" unless (@ARGV == 3 || @ARGV == 4);

my ($outdir, $prefix) = @ARGV[1,2];


my $count = 1;
my $mutProSeq = "";
my $singleProLen = 500;
my $gap = "BBBBBBBBBB";

open (OUT_FA, "> $ARGV[1]") || die $!;
if (@ARGV == 4) {
    ($ARGV[3] =~ /\.gz$/)?(open (IN, "gzip -cd $ARGV[3] | ") || die $!):(open (IN, "$ARGV[3]") || die $!);
    while (<IN>) {
        chomp;
        print OUT_FA $_."\n";
    }
    close IN;
}

open (OUT_POS, "> $ARGV[2]") || die $!;
print OUT_POS join("\t", qw(MUTPROid  startInPro  endInPro  MUTPEPid  MUTPEPseq))."\n";
open (IN, "$ARGV[0]") || die $!;
$/ = ">"; <IN>;
while (<IN>) {
    chomp;
    my @inf = split /\n/, $_;
    my $pepID = $inf[0];
    my $pepSeq = $inf[1];
    if (!$mutProSeq) {
        $mutProSeq = $pepSeq;
    }else {
        $mutProSeq .= "$gap$pepSeq";
    }
    my $start = length($mutProSeq)-length($pepSeq)+1;
    my $end = length($mutProSeq);
    my $mutProID = sprintf("MUTPRO%06d", $count);
    print OUT_POS join("\t", $mutProID, $start, $end, $pepID, $pepSeq)."\n";
    if (length($mutProSeq) > $singleProLen) {
        print OUT_FA ">$mutProID\n$mutProSeq\n";
        $mutProSeq = "";
        $count ++;
    }
}
close IN;
close OUT_POS;
close OUT_FA;
