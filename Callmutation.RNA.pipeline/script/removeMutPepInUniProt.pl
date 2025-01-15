#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <mutant.peptides.raw.gz>   <UniProtDB.fa[.gz]>   <outdir>   <prefix>\n\n" unless (@ARGV == 4);

my ($outDir, $outPrefix) = @ARGV[2,3];

($ARGV[0] =~ /\.gz$/)?(open (MUT, "gzip -cd $ARGV[0] | ") || die $!):(open (MUT, "$ARGV[0]") || die $!);
my %head = ();
my @tmp = split /\t/, <MUT>; chomp $tmp[-1];
for (my $i=0; $i<@tmp; $i++) {
    $head{$tmp[$i]} = $i;
}

my %mutAA;
while (<MUT>) {
    chomp;
    my @inf = split /\t/, $_;
    my $pep = $inf[$head{"mutAA"}];
    my $pepHead = substr($pep, 0, 4);
    $pepHead =~ s/I/L/g;
    $mutAA{$pepHead}{$pep} = 1;
}
close MUT;

($ARGV[1] =~ /\.gz$/)?(open (DBF, "gzip -cd $ARGV[1] | ") || die $!):(open (DBF, "$ARGV[1]") || die $!);
$/ = ">"; <DBF>;

my %pep2pro;
my ($block, $num) = (0, 0);
while (<DBF>) {
    chomp;
    my @inf = split /\n/, $_;
    my $id = (split /\s+/, $inf[0])[0];
    my $seq = join("", @inf[1..$#inf]);
    $seq =~ s/I/L/g;
    foreach my $pepHead (sort{$a cmp $b} keys %mutAA) {
        next if ($seq !~ /$pepHead/);
        foreach my $pep (sort{$a cmp $b} keys %{$mutAA{$pepHead}}) {
            my $transPEP = $pep;
            $transPEP =~ s/I/L/g;
            next if ($seq !~ /$transPEP/);
            push @{$pep2pro{$pep}}, $id;
        }
    }
    $num ++;
    $block ++;
    if ($block == 500) {
        my $date = `date`; chomp $date;
        print "$id\t$num\t$date";
        $block = 0;
    }
}
close DBF;
$/ = "\n";

open (OUT, " | gzip > $outDir/$outPrefix.SnvIndel.mutant.peptides.noDBmatch.gz") || die $!;
print OUT join("\t", qw(Chr  Start  End  Ref  Alt  Region  Gene  GeneDetail  mutType  mutPEPlen  mutDetail  mutStrand  wildAA  mutAA  wildNT  mutNT))."\n";

($ARGV[0] =~ /\.gz$/)?(open (MUT, "gzip -cd $ARGV[0] | ") || die $!):(open (MUT, "$ARGV[0]") || die $!);
%head = ();
@tmp = split /\t/, <MUT>; chomp $tmp[-1];
for (my $i=0; $i<@tmp; $i++) {
    $head{$tmp[$i]} = $i;
}

my (%mutSite, %regionPepStat, %typePepStat, %uniquePepStat, %filterPepStat);
while (<MUT>) {
    chomp;
    my @inf = split /\t/, $_;
    my $pep = $inf[$head{"mutAA"}];
    #next if (exists $pep2pro{$pep});
    if (exists $pep2pro{$pep}) {
        $filterPepStat{$pep} = 1;
        #print join("\t", @inf, join(",", @{$pep2pro{$pep}}))."\n";
        next;
    }
    print OUT join("\t", @inf)."\n";
    
    my ($chr, $start, $end, $region, $type) = @inf[$head{"Chr"}, $head{"Start"}, $head{"End"}, $head{"Region"}, $head{"mutType"}];
    my $site = join("-", $chr, $start, $end);
    $mutSite{$site} = 1;
    $regionPepStat{$region}{$pep} = 1;
    $typePepStat{$type}{$pep} = 1;
    $uniquePepStat{$pep} = 1;
}
close MUT;

open (STAT, "> $outDir/$outPrefix.SnvIndel.noDBmatch.mutPepCount.stat.tsv") || die $!;
my $totalUniqPep = scalar(keys %uniquePepStat);
my $mutSiteCounts = scalar(keys %mutSite);
my $ratio = sprintf("%.2f", $totalUniqPep/$mutSiteCounts);
print STAT "#Number of total mutant sites from $outPrefix: $mutSiteCounts\n";
print STAT "#Number of total unique peptides from $outPrefix: $totalUniqPep\n";
print STAT "#Number of peptides generated per mutant site: $ratio\n\n";

my $filterNum = scalar(keys %filterPepStat);
print STAT "Number of filter peptides with DB match: $filterNum\n";
print STAT "#Summary of peptide number by gene region\n";
foreach my $key (sort{$a cmp $b} keys %regionPepStat) {
    my $regionUniq = scalar(keys %{$regionPepStat{$key}});
    my $rate = sprintf("%.4f", $regionUniq/$totalUniqPep*100);
    print STAT join("\t", $key, $regionUniq, $rate)."\n";
}
print STAT "\n#Summary of peptide number by mutant type\n";
foreach my $key (sort{$a cmp $b} keys %typePepStat) {
    my $typeUniq = scalar(keys %{$typePepStat{$key}});
    my $rate = sprintf("%.4f", $typeUniq/$totalUniqPep*100);
    print STAT join("\t", $key, $typeUniq, $rate)."\n";
}
close STAT;
