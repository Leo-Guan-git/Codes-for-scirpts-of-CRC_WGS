#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);

=head1 Description

=head1 Version

        Version: 0.1    Date: 2022-04-02   Author: xianghaitao@genomics.cn

=head1 Options

        -input             <s> : path to input file to split (mandatory)
        -outdir            <s> : path to output directory [default: pwd]
        -prefix            <s> : path to output directory [default: fileSplit]
        -header            <s> : whether there is a header in the input file, yes/no [default: yes]
        -inputFileLines    <n> : the number of lines of the input file [default: '0' means not specify]
        -splitFileNum      <n> : the number of subfiles the input file will be split into [default: 20]
        -singleFileLines   <n> : the number of lines per subfile into which the input file will be split [default: '0' means not specify]
        -help                  : show this help information and exit.

=cut

my ($input, $outdir, $prefix, $header, $inputFileLines, $splitFileNum, $singleFileLines, $help);
GetOptions (
    "input:s" => \$input,
    "outdir:s" => \$outdir,
    "prefix:s" => \$prefix,
    "header:s" => \$header,
    "inputFileLines:n" => \$inputFileLines,
    "splitFileNum:n" => \$splitFileNum,
    "singleFileLines:n" => \$singleFileLines,
    "help" => \$help
);

die `pod2text $0` if (!defined $input || $help);

$outdir ||= `pwd`;
chomp $outdir;
$outdir = abs_path($outdir);
system("mkdir -p $outdir");

$prefix ||= "fileSplit";
$header ||= "yes";

$splitFileNum ||= 20;
$inputFileLines ||= 0;
$singleFileLines ||= 0;

if ($singleFileLines == 0) {
    if ($inputFileLines == 0) {
        if ($input =~ /.gz$/) {
            $inputFileLines = `gzip -cd $input | wc -l`;
            chomp $inputFileLines;
        }else {
            $inputFileLines = `wc -l $input`;
            chomp $inputFileLines;
            $inputFileLines = (split /\s+/, $inputFileLines)[0];
        }
        $inputFileLines -= 1 if ($header eq 'yes');
    }
    $singleFileLines = int($inputFileLines/$splitFileNum) + 1;
}

my $base = basename($input);
print "### File cutting statistics for $base:\n";
print join("\t", qw(inputFileLines  splitFileNum  singleFileLines))."\n";
print join("\t", $inputFileLines,  $splitFileNum, $singleFileLines)."\n";

($input =~ /\.gz$/)?(open (IN, "gzip -cd $input | ") || die $!):(open (IN, "$input") || die $!);
my ($line, $block) = (0, 1);
my $subFileID = sprintf("%04d", $block);
open (OUT, "> $outdir/$prefix.$subFileID.txt") || die $!;
my @head = ();
if ($header =~  /yes/i) {
    @head = split /\t/, <IN>; chomp $head[-1];
    print OUT join("\t", @head)."\n";
}
while (<IN>) {
    chomp;
    $line ++;
    print OUT $_."\n";
    if ($line == $singleFileLines) {
        close OUT;
        
        $block ++;
        $subFileID = sprintf("%04d", $block);
        open (OUT, "> $outdir/$prefix.$subFileID.txt") || die $!;
        print OUT join("\t", @head)."\n" if ($header =~  /yes/i);
        $line = 0;
    }
}
close IN;
close OUT;
