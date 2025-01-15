#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <SnvIndel.mutant.peptides.noDBmatch.withORFmatch.unique[.gz]>   <outdir>   <prefix>\n\n" unless (@ARGV == 3);

my ($outdir, $prefix) = @ARGV [1,2];

($ARGV[0] =~ /\.gz$/)?(open(IN, "gzip -cd $ARGV[0] | ") || die $!):(open(IN, "$ARGV[0]") || die $!);
<IN>;
my %hash;
while (<IN>) {
    chomp;
    my @inf = split /\t/, $_;
    push @{$hash{$inf[1]}}, $inf[0];
}
close IN;

open (STAT, "> $outdir/$prefix.peptide.counts.stat") || die $!;
print STAT join("\t", qw(setAlen  setBlen  setAnum  setBnum))."\n";
my @len = sort{$a <=> $b} keys %hash;
for (my $i=0; $i<@len; $i++) {
    
    if ($i+1 < @len) {
        my @setB = ();
        for (my $j=$i+1; $j<@len; $j++) {
            push @setB, @{$hash{$len[$j]}};
        }
        
        my $setAnum = scalar(@{$hash{$len[$i]}});
        my $setBnum = scalar(@setB);
        
        my $setBlen = join("-", @len[$i+1, $#len]);
        $setBlen = $len[-1] if ($i+2 == @len);
        
        my $subdir = "$outdir/length_$len[$i].vs.$setBlen";
        `mkdir -p $subdir`;
        open (OUT1, "> $subdir/$prefix.length$len[$i].pep.list") || die $!;
        print OUT1 join("\n", @{$hash{$len[$i]}})."\n";
        close OUT1;
        
        open (OUT2, "> $subdir/$prefix.length$setBlen.pep.list") || die $!;
        print OUT2 join("\n", @setB)."\n";
        close OUT2;
        
        print STAT join("\t", $len[$i], $setBlen, $setAnum, $setBnum)."\n";
    }else {
        my $subdir="$outdir/length_$len[$i]";
        `mkdir -p $subdir`;
        open (OUT1, "> $subdir/$prefix.length$len[$i].retained.peptides.list") || die $!;
        print OUT1 join("\n", @{$hash{$len[$i]}})."\n";
        close OUT1;
    }
}
close STAT;
