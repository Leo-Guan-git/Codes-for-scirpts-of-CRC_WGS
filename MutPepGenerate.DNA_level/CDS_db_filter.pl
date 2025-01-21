#! /usr/bin/env perl -w
# @Author: Xiangyu Guan
# @Date:   2025-01-21
# @Last Modified by:   Xiangyu Guan
# @Last Modified time: 2025-01-21
use strict;
use warnings;
use POSIX;
use File::Basename;
use Getopt::Std;
use vars qw($opt_r $opt_i $opt_p);
getopts('r:i:p:');

my ($ref_f, $in_f, $out_prefix, $out_f);
my %ref_seqs;
$ref_f = "/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/MS/MQ.reference.fa";

if(! ($opt_i || $opt_p)){
    print "Usage: $0\n";
    print "--------------------------------------------------------------------\n";
    print "-r uniprot protein databases.............[Homo.sapiens.9606.fa[.gz]]\n";
    print "-i input fasta file......................<*.vcf[.gz]>\n";
    print "-p output_prefix.........................<./test>\n";
    die "Please Check your parameters\n";    
}

$ref_f = $opt_r if($opt_r);
$in_f = $opt_i if($opt_i);
$out_prefix = $opt_p if($opt_p);
$out_f = $out_prefix.".protein.DB.fasta.gz";

Message("Reference file: $ref_f");
Message("Input file: $in_f");
Message("Output file: $out_f");

Message("Reference Parsing...");
if($ref_f =~ /\.gz$/){ open REF, "gzip -dc $ref_f |" or die $!; }else{ open REF, "<$ref_f" or die $!; }
open OUT, "| gzip -c > $out_f" or die $!;
my ($line, $head, $seq, $db, $uniq_id, $Entry_name, $new_line);
my %pep;
while ($line = <REF>){
    chomp $line;
    if ($line =~ /^>/){
        if($seq){
            $ref_seqs{$uniq_id} = $seq;
            if ($uniq_id !~ /Biognosys/){ $pep{Origin} += 1; }
            print OUT "$new_line\n";
            print OUT SeqLineSplit($seq)."\n";
        }
        $head = (split /\s+/, $line)[0];
        ($db, $uniq_id, $Entry_name) = split /\|/, $head;
        $new_line = $line;
        $seq = "";
    }else{
        $seq .= $line;
    }
}
if($seq){
    $ref_seqs{$uniq_id} = $seq;
    if ($uniq_id !~ /Biognosys/){ $pep{Origin} += 1; }
    print OUT "$new_line\n";
    print OUT SeqLineSplit($seq)."\n";
}
close REF;
Message("Reference Proteins Parsing done");
Message("$pep{Origin} reference proteins are selected");

Message("Mutated Proteins filtering...");
if($in_f =~ /\.gz$/){ open IN, "gzip -dc $in_f |" or die $!; }else{ open IN, "<$in_f" or die $!; }
$seq = "";
while ($line = <IN>) {
    chomp $line;
    if ($line =~ /^>/){
        if($seq){
            my $included = 0;
            foreach my $id(keys %ref_seqs){
                if ($ref_seqs{$id} =~ /$seq/){ $included = 1; last; }
            }
            unless($included){
                $pep{Mutated} += 1;
                print OUT $new_line;
                print OUT SeqLineSplit($seq)."\n";
            }
        }
        $line =~ /^>mut\|([\w\.]+)_prot_([\w\.]+)_(mut\d+) /;
        my ($chrom, $pro_id, $mut_id) = ($1, $2, $3);
        $new_line = ">mut|".$mut_id."|".$chrom."_prot_".$pro_id." ".$'."\n";
        $seq = "";
    }else{
        $seq .= $line;
    }
}
if($seq){
    my $included = 0;
    foreach my $id(keys %ref_seqs){
        if ($ref_seqs{$id} =~ /$seq/){ $included = 1; last; }
    }
    unless($included){
        $pep{Mutated} += 1;
        print OUT $new_line;
        print OUT SeqLineSplit($seq)."\n";
    }
}

close IN;
close OUT;
Message("Mutated Proteins filtering done");
Message("$pep{Mutated} mutated proteins are selected");
Message("Jons are done");
sub SeqLineSplit{
    my $seq = $_[0];
    my $LINEBASES = 80;
    my @seq_splited;
    for (my $n = 0; $n < length($seq); $n += $LINEBASES) {
        push @seq_splited, substr($seq, $n, $LINEBASES);
    }
    return join("\n", @seq_splited);
}

sub Message{
    my $message = $_[0];
    my $datestring = localtime();
    $message = "$datestring >>> $message\n";
    print $message;
}