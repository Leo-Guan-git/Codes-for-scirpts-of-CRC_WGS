#! /usr/bin/env perl -w
# @Author: Xiangyu Guan
# @Date:   2025-01-21
# @Last Modified by:   Xiangyu Guan
# @Last Modified time: 2025-01-21
use strict;
use warnings;
use POSIX;
use File::Basename;
use POSIX;
use Fcntl qw( SEEK_SET );

@ARGV ==5 or die "Usage: perl $0 <orfipy.fa[.gz]> <MutPep.tsv> <Outpath> <reference.fa[.gz]> <max_pep_length>\n";
my ($orffile, $MutPep, $outfile, $reference, $max_pep_length) = @ARGV;
my ($chrom, $orfid, $ORF_start, $ORF_stop, $STRAND, $length, $CONDON_start, $CONDON_stop);
my ($orf_num, $FRAME, $bin_id);
my (%OrfIntervals, %bin_size, %pep_num, %selected_pep);

sub Message{
    my $message = $_[0];
    my $datestring = localtime();
    print "\t $datestring >>> $message\n";
}

Message("Mut-Peptides filter Started");
Message("ORF predicted file: $orffile");
Message("Mut-Peptide file: $MutPep");
Message("Output file: $outfile");
Message("Reference file: $reference");
Message("Max Peptide length: $max_pep_length");

# referece fai file Parse:
my $reference_fai = $reference.".fai";
my %faidx;
if (! -e $reference_fai){
    # generate fai file if not exist
    system("samtools faidx $reference")
}

open FAI, "<$reference_fai";
while (<FAI>){
    chomp $_;
    my ($CHROM, $LENGTH, $OFFSET, $LINEBASES, $LINEWIDTH) = split /\t/, $_;
    $faidx{$CHROM}{"LENGTH"} = $LENGTH; 
    $faidx{$CHROM}{"OFFSET"} = $OFFSET; 
    $faidx{$CHROM}{"LINEBASES"} = $LINEBASES; 
    $faidx{$CHROM}{"LINEWIDTH"} = $LINEWIDTH;
}
close FAI;

# in order to use seek function in read-in target sequence,
# .gz file must be decompressed.
my $tmp_file;
if ($reference =~ /(\S+)\.gz$/){
    $tmp_file = $1;
    system("gzip -dc $reference > $tmp_file");
    open REF, "< $tmp_file";
}else{
    open REF, "< $reference";
}

sub FastaSeqRead{
    # inorder to fit for complete human genome.
    # do not load a big reference to memory.
    my $CHROM = $_[0];
    my $seq;
    my $start_byte = $faidx{$CHROM}{"OFFSET"};
    my $length_byte = floor($faidx{$CHROM}{"LENGTH"} / $faidx{$CHROM}{"LINEBASES"}) + $faidx{$CHROM}{"LENGTH"};
    seek REF, $start_byte, SEEK_SET or die $!;
    read REF, $seq, $length_byte;
    $seq =~ tr/\n//d;
    $seq =~ tr/a-z/A-Z/;
    return $seq
}

Message("chromosome segementation start");

# bin size calculate
# for orfipy predicted orf result, matching steps can be
# sharply reduced by split to 6 different frame parts
# as well as segement each frame based on maximum length
# of ORF regions.
$orf_num = 0;
if ($orffile =~ /\.gz$/){open ORF, "gzip -dc $orffile |" or die $!; }else{ open ORF, "< $orffile" or die $!; }
while (my $line=<ORF>){
    if ($line !~ /^\>/){ next; }
    $orf_num += 1;
    chomp $line;
    $line =~ /\>(chr\S+)_(ORF\.\d+) \[(\d+)\-(\d+)\]\(([+-])\) \S+ length\:(\d+) frame\:\-?[123] start\:([ATCG]{3}) stop\:([ATCG]{3})/;
    ($chrom, $orfid, $ORF_start, $ORF_stop, $STRAND, $length, $CONDON_start, $CONDON_stop) = @{^CAPTURE};
    $FRAME = $ORF_start % 3;
    if (!exists($bin_size{$chrom}{$STRAND.$FRAME})){
        $bin_size{$chrom}{$STRAND.$FRAME} = $length+1;
    }else{
        if ($length>$bin_size{$chrom}{$STRAND.$FRAME}){
            $bin_size{$chrom}{$STRAND.$FRAME}=$length+1;
        }
    }
}
close ORF;

Message("chromosome segementation complete");

for my $chrom(keys %bin_size){
    for my $frame(keys %{$bin_size{$chrom}}){
        if ($bin_size{$chrom}{$frame} > 100000){
            # for ORF prediction result from orfipy, ORF regions could be too large 
            # to reduce the segementation effects.
            $bin_size{$chrom}{$frame} = 100000;
        }
        Message("$chrom Frame $frame bin size: $bin_size{$chrom}{$frame}")
    }
}

Message("ORF Parsing start");

# ORF message Parsing
if ($orffile =~ /\.gz$/){open ORF, "gzip -dc $orffile |" or die $!; }else{open ORF, "< $orffile" or die $!;}
$orf_num=0;
while (my $line=<ORF>){
    if ($line !~ /^\>/){next;}
    $orf_num += 1;
    chomp $line;
    $line =~ /\>(chr\S+)_(ORF\.\d+) \[(\d+)\-(\d+)\]\(([+-])\) \S+ length\:(\d+) frame\:\-?[123] start\:([ATCG]{3}) stop\:([ATCG]{3})/;
    ($chrom, $orfid, $ORF_start, $ORF_stop, $STRAND, $length, $CONDON_start, $CONDON_stop) = @{^CAPTURE};
    $ORF_start += 1;
    $FRAME = $ORF_start % 3;
    $bin_id = ceil( $ORF_start / $bin_size{$chrom}{$STRAND.$FRAME});
    $OrfIntervals{$chrom}{$STRAND.$FRAME}{$bin_id}{"$ORF_start\-$ORF_stop"} = join("\t", ($orfid, $ORF_start, $ORF_stop, $CONDON_start, $CONDON_stop));
}
close ORF;

Message("ORF Parsing Done");
Message("ORF number: $orf_num");
Message("Mutpeptide filtering start");

if ($MutPep =~ /\.gz$/){ open PEP, "gzip -dc $MutPep |" or die $!; }else{ open PEP, "< $MutPep"; }
if ($outfile =~ /\.gz$/){open TSV, "| gzip > $outfile" or die $!; }else{open TSV, "> $outfile" or die $!; }
my $header = readline PEP;
my @header = split /\t/, $header;
my @pep_length;
for my $index(@header){
    if ($index =~ /^PEP_(\d+)$/){
        push @pep_length, $1;
    }
}

print TSV "CHROM\tPOS\tREF\tNEW_POS\tALT\tis_somatic\tMUT_TYPE\tSTRAND\tNS_START\tPEP\tNUCLE\tPEP_ID\tORF_start\tORF_stop\tCONDON_START\tCONDON_END\n";
my $processing_pos = 0;
# my %selected_pep;
my $max_ns_length = $max_pep_length * 3;
while (my $line=<PEP>){
    chomp $line;
    my ($CHROM, $POS, $REF, $NEW_POS, $ALT, $is_somatic, $MUT_TYPE, @Peps_info) = split /\t/, $line;
    my $ref_seq = FastaSeqRead($CHROM);
    my %Peps;
    my $NEW_END = $NEW_POS + length($ALT) - 1;
    for my $n(0..$#Peps_info){
        $length = $pep_length[$n];
        $Peps_info[$n] =~ /\(\[([\s\S]+)\]\, \[([\s\S]+)\]\, \[([\s\S]+)\]\)/;
        my ($strand, $start, $pep) = ($1,$2,$3);
        $strand =~ tr/\'//d;
        $pep =~ tr/\'//d;
        my @strand = split /\, /, $strand;
        my @start = split /\, /, $start;
        my @pep = split /\, /, $pep;
        for my $i(0..$#strand){
            if ($strand[$i] =~ /\s/){next;}
            my $STRAND = $strand[$i];
            my $NS_START = $start[$i];
            my $PEP = $pep[$i];
            my $NUCLE = substr $ref_seq, $NS_START-1, $length*3;
            my $FRAME = $NS_START % 3;
            my $NS_END = $NS_START + $length*3 - 1;
            $Peps{$CHROM}{"$NEW_POS\-$NEW_END"}{$STRAND.$FRAME}{"$NS_START\-$NS_END"} = join("\t", ($CHROM, $POS, $REF, $NEW_POS, $ALT, $is_somatic, $MUT_TYPE, $STRAND, $NS_START, $PEP, $NUCLE));
            $pep_num{$length} += 1;
        }
    }
    $processing_pos += 1;
    if (($processing_pos % 10000)==0){Message("$processing_pos mutate positions are processed");}
    for my $CHROM(keys %Peps){
        for my $pos(sort keys %{$Peps{$CHROM}}){
            my ($NEW_POS, $NEW_END) = split /\-/, $pos;
            # for my $NEW_END(sort {$a <=> $b} keys %{$Peps{$CHROM}{$NEW_POS}}){
            for my $frame(sort keys %{$Peps{$CHROM}{$pos}}){
                my $found = 0;
                my $bin_id = ceil($NEW_POS / $bin_size{$CHROM}{$frame});
                if (exists($OrfIntervals{$CHROM}{$frame}{$bin_id})){
                    for my $ORF_pos(sort keys %{$OrfIntervals{$CHROM}{$frame}{$bin_id}}){
                        if ($found == 1){last;}
                        my ($ORF_start, $ORF_stop) = split /\-/, $ORF_pos;
                        # for my $ORF_stop(sort keys %{$OrfIntervals{$CHROM}{$frame}{$bin_id}{$ORF_start}}){
                        if ($found == 1){last;}
                        if (($ORF_start <= $NEW_POS) && ($NEW_END <= $ORF_stop)){
                            if ((($NEW_POS - $ORF_start + 1) <= $max_ns_length) && (($ORF_stop - $NEW_END + 1) >= $max_ns_length)){
                                for my $NS_pos(sort keys %{$Peps{$CHROM}{$pos}{$frame}}){
                                    my ($NS_START, $NS_END) = split /\-/, $NS_pos;
                                    # for my $NS_END(sort {$a <=> $b} keys %{$Peps{$CHROM}{$pos}{$frame}{$NS_START}}){
                                    print TSV "$Peps{$CHROM}{$pos}{$frame}{$NS_pos}\t$OrfIntervals{$CHROM}{$frame}{$bin_id}{$ORF_pos}\n";
                                    $found = 1;
                                    $selected_pep{$NS_END-$NS_START+1} += 1;
                                }
                            }else{
                                for my $NS_pos(sort keys %{$Peps{$CHROM}{$pos}{$frame}}){
                                    my ($NS_START, $NS_END) = split /\-/, $NS_pos;
                                    # for my $NS_END(sort {$a <=> $b} keys %{$Peps{$CHROM}{$pos}{$frame}{$NS_START}}){
                                    if (($ORF_start <= $NS_START) && ($NS_END <= $ORF_stop)){
                                        print TSV "$Peps{$CHROM}{$pos}{$frame}{$NS_pos}\t$OrfIntervals{$CHROM}{$frame}{$bin_id}{$ORF_pos}\n";
                                        $found = 1;
                                        $selected_pep{$NS_END-$NS_START+1} += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                $bin_id = $bin_id - 1;
                if (($found==0) && exists($OrfIntervals{$CHROM}{$frame}{$bin_id})){
                    for my $ORF_pos(sort keys %{$OrfIntervals{$CHROM}{$frame}{$bin_id}}){
                        if ($found == 1){last;}
                        my ($ORF_start, $ORF_stop) = split /\-/, $ORF_pos;
                        # for my $ORF_stop(sort {$a <=> $b} keys %{$OrfIntervals{$CHROM}{$frame}{$bin_id}{$ORF_start}}){
                        if ($found == 1){last;}
                        if (($ORF_start <= $NEW_POS) && ($NEW_END <= $ORF_stop)){
                            if ((($NEW_POS - $ORF_start + 1) <= $max_ns_length) && (($ORF_stop - $NEW_END + 1) >= $max_ns_length)){
                                for my $NS_pos(sort keys %{$Peps{$CHROM}{$pos}{$frame}}){
                                    my ($NS_START, $NS_END) = split /\-/, $NS_pos;
                                    # for my $NS_END(sort {$a <=> $b} keys %{$Peps{$CHROM}{$pos}{$frame}{$NS_START}}){
                                    print TSV "$Peps{$CHROM}{$pos}{$frame}{$NS_pos}\t$OrfIntervals{$CHROM}{$frame}{$bin_id}{$ORF_pos}\n";
                                    $found = 1;
                                    $selected_pep{$NS_END-$NS_START+1} += 1;
                                }
                            }else{
                                for my $NS_pos(sort keys %{$Peps{$CHROM}{$pos}{$frame}}){
                                    my ($NS_START, $NS_END) = split /\-/, $NS_pos;
                                    # for my $NS_END(sort {$a <=> $b} keys %{$Peps{$CHROM}{$pos}{$frame}{$NS_START}}){
                                    if (($ORF_start <= $NS_START) && ($NS_END <= $ORF_stop)){
                                        print TSV "$Peps{$CHROM}{$pos}{$frame}{$NS_pos}\t$OrfIntervals{$CHROM}{$frame}{$bin_id}{$ORF_pos}\n";
                                        $found = 1;
                                        $selected_pep{$NS_END-$NS_START+1} += 1;
                                    }
                                }
                            }
                        }
                        # }
                    }
                }
            }
        }
    }
}
close REF;
close PEP;
close TSV;

Message("Mutate positions filtering are done");
if (-e $tmp_file){system("rm $tmp_file");}
for my $length(sort {$a <=> $b} keys %pep_num){
    my $ns_length = $length * 3;
    my $percentage = $selected_pep{$ns_length} / $pep_num{$length} * 100;
    $percentage = sprintf "%.3f", $percentage;
    Message("$selected_pep{$ns_length} of $pep_num{$length} ( $percentage\% ) raw peptides for $length long peptides");
}

Message("Mission Complete!")