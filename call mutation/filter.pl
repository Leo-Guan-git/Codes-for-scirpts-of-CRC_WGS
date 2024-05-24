#!usr/bin/perl -w
use strict;

if(@ARGV<5)
{
	print "$0 snv.softmerge.raw.snv.xls  snv.softmerge.raw.indel.xls mutect2.filter2.vcf.gz mutect.raw.txt.gz snv.softmerge.filter.snv.xls snv.softmerge.filter.indel.xls\n";
	exit 1;
}

my ($snvf,$indelf,$mutect2,$mutect,$snvout,$indelout)=($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5]);
my %mutect2=();
my %mutect=();
open IN,"gzip -dc $mutect2|" or die $!;
while(<IN>)
{
	chomp;
	next if(/^#/);
	my @a=split /\t/,$_;
	if($a[6] ne 'PASS')	
	{
		my @alts=split /\,/,$a[4];
		foreach my $alt(@alts)
		{
			my $key="$a[0]\t$a[1]\t$a[3]\t$alt";
			$mutect2{$key}=1;
		}
	}
}
close IN;

open IN,"gzip -dc $mutect|" or die $!;
while(<IN>)
{
	chomp;
	next if(/^#/);
	my @a=split /\t/,$_;
	if($a[50] ne 'KEEP')	
	{
		my $key="$a[0]\t$a[1]\t$a[3]\t$a[4]";
		$mutect{$key}=1;
	}
}
close IN;

open IN,$snvf or die $!;
open OUT,">$snvout" or die $!;
while(<IN>)
{
	chomp;
	my @a=split /\t/,$_;
	my $line=$_;
	my $key="$a[0]\t$a[1]\t$a[3]\t$a[4]";
	next if(exists $mutect2{$key} || exists $mutect{$key});
	print OUT $line."\n";
}
close IN;
close OUT;

open IN,"$indelf" or die $!;
open OUT,">$indelout" or die $!;
while(<IN>)
{
	chomp;
	my @a=split /\t/,$_;
	my $line=$_;
	my $key="$a[0]\t$a[1]\t$a[3]\t$a[4]";
	next if(exists $mutect2{$key});
	print OUT $line."\n";
}
close IN;
close OUT;

