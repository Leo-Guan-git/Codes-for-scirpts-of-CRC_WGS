#!/usr/bin/perl -w
# @Author: guanxiangyu
# @Date:   2023-02-18
# @Last Modified by:   guanxiangyu
# @Last Modified time: 2023-11-17 16:17:47
use strict;
use POSIX;
use warnings;
#use List::MoreUtils qw( minmax );
use File::Basename;
use Data::Dumper;
######################### Update Info #################################
# 2023.02.23 by guanxiangyu@genomics.cn
# 1. support PE fastq or bam input, SE fastq can not by recorected
# 2. use STAR for alignment
# 3. add readgroup contents if it is missed in input bam files 
# 4. use GATK for BAM file post-process 
# 5. tumor only mode are imcomplete
######################## Command paraters getting #####################
my $script = basename $0;
@ARGV == 3 or die "Usage: perl $script <input.list> <Configs> <Outdir>\n";
my ($inputlist,$config,$outdir) = @ARGV;
system("mkdir -p $outdir") unless (-e "$outdir");
#######################################################################
my (%configs,%scripts);
my (%fqinputs,%cleanfqs,%splitfqs,%cleanfastq,%correctfastq,%phred,%hlatype,%unmapfqs);
my (%baminputs,%realnbam,%realnbam2,%bqsrbam,%rnabam,%SplitedBams);
my (%snvvcf,%indelvcf,%fusionresult,%exprresult);

########################### Main Parameters ###########################

readConfigs($config);

readList($inputlist);
$configs{bampostprocess} = "true" if($configs{alignment} =~ /true/i);
$configs{addreadgroup} = "false" if ($configs{alignment} =~ /true/i);
if($configs{tumor_only} =~ /true/i){
    $configs{somatic} = "false";
    $configs{svaba} = "false";
    $configs{MuSE} = "false";
    $configs{strelka} = "false";
    $configs{strelka2} = "false";
    $configs{mutect} = "false";
    $configs{mutect2} = "false";
    $configs{bamqc} = "false";
}
my @chrs;
open LIST,"<$configs{reference_block}" or die "$!";
while(my $line = <LIST>){
    chomp $line;
    push @chrs, $line;
}
close LIST;
makeDir(%fqinputs) if($configs{alignment} =~ /true/i);
makeDir(%baminputs) if($configs{alignment} =~ /false/i);
runsplit(%fqinputs) if ($configs{alignment} =~ /true/i and $configs{split_num} > 1);

if ($configs{fastqclean} =~ /true/i){
    if($configs{split_num} > 1){ FastqClean(%splitfqs); }else{ FastqClean(%fqinputs); }
}else{
    FastqLink(%fqinputs);
}
FastqCorrect(%cleanfastq) if($configs{fastqcorrect} =~ /true/i);
##align BQSR HaplotypeCaller-Sentieon
#rmdup + BQSR + haplotyper - Sentieon
if($configs{alignment} =~ /true/i){
    if($configs{fastqcorrect} =~ /true/i){ Alignment(%correctfastq); }else{ Alignment(%cleanfastq); }
}
else{
    if($configs{addreadgroup} =~ /true/i){
        Alignment(%baminputs);
    }
}
##QC
bamqc(%bqsrbam) if($configs{bamqc} =~ /true/i);

##call SV
svaba(%bqsrbam) if($configs{svaba} =~ /true/i);
GATK_TumorOnly(%bqsrbam) if($configs{tumor_only} =~ /true/i);

## somatic mutation calling
mutectSNV(%bqsrbam) if ($configs{mutect} =~ /true/i);
mutect2(%bqsrbam) if ($configs{mutect2} =~ /true/i);
strelka2(%bqsrbam) if ($configs{strelka2} =~ /true/i);
muse(%bqsrbam) if ($configs{MuSE} =~ /true/i);
strelka(%bqsrbam) if ($configs{strelka} =~ /true/i);
mergeSnvInDel(%bqsrbam) if ($configs{mutect} =~ /true/i and $configs{mutect2} =~ /true/i);

#HLA type
#HLAtyping()
if(-e $configs{hlalist}){
    readHLAlist();
}else{
    if($configs{hlasomatic} =~ /true/i){
        hlasomatic(%bqsrbam);
    }
}

# Kraken2
kraken2(%unmapfqs) if($configs{kraken} =~ /true/i);

# Generate day info
chomp(my $day = `date +%Y%m%d`);
## run scripts
edgeList($day);
## Generate run.sh
open RUNSH, ">$outdir/shell_run/run.$configs{projectname}.$day.sh" or die "$!";
print RUNSH "$configs{monitor} taskmonitor --q1 $configs{queue} --q2 $configs{queue2}";
print RUNSH " --P1 $configs{priority} --P2 $configs{priority2}";
print RUNSH " -p $configs{projectname} -i $outdir/shell_run/edge.$configs{projectname}.$day.list -f 1\n";
close RUNSH;


########################### Functions ###########################
sub readConfigs{
    my ($configfile) = @_;
    open CFG, "<$configfile" or die "$!";
    while(my $line = <CFG>){
        chomp $line;
        next if ($line =~ /^#/);
        $configs{$1} = $2 if ($line =~ /^(.*)=(.*);/);
    }
    close CFG;
}

sub readList{
    my ($inputs) = @_;
    open INPUT, "<$inputs" or die "$!";
    my $num = 0;
    while(my $line = <INPUT>){
        chomp $line;
        my ($sampid,$dtype,$phred,$lane,$input) = split /\s+/,$line;
        $dtype = "Tumor" if($dtype =~ /Tumor/i);
        $dtype = "Normal" if($dtype =~ /Normal/i);
        if ($configs{alignment} =~ /false/i){
            ## bam file input
            push @{$baminputs{$sampid}{$dtype}{bam}}, $input; 
        }else{
            # fastq file input
            my ($fq1,$fq2) = split /,/, $input;
            push @{$fqinputs{$sampid}{$dtype}{$lane}{fq1}}, $fq1;
            push @{$fqinputs{$sampid}{$dtype}{$lane}{fq2}}, $fq2;
        }
        $phred{$sampid}{$dtype} = $phred;
    }
    close INPUT;
}

sub makeDir{
    my (%hash) = @_;
    system("mkdir -p $outdir/shell_run") unless (-e "$outdir/shell_run");
    for my $sampid(sort keys %hash){
        system("mkdir -p $outdir/$sampid/0.shell") unless (-e "$outdir/$sampid/0.shell");
        system("mkdir -p $outdir/$sampid/1.clean") unless (-e "$outdir/$sampid/1.clean");
        system("mkdir -p $outdir/$sampid/2.alignment") unless (-e "$outdir/$sampid/2.alignment");
        system("mkdir -p $outdir/$sampid/3.somatic") unless (-e "$outdir/$sampid/3.somatic");
        system("mkdir -p $outdir/$sampid/4.microbiome") unless (-e "$outdir/$sampid/4.microbiome");
        ## shell
        system("mkdir -p $outdir/$sampid/0.shell/a01.clean") unless (-e "$outdir/$sampid/0.shell/a01.clean");
        system("mkdir -p $outdir/$sampid/0.shell/a02.alignment") unless (-e "$outdir/$sampid/0.shell/a02.alignment");
        system("mkdir -p $outdir/$sampid/0.shell/a04.mutect") unless (-e "$outdir/$sampid/0.shell/a04.mutect");
        system("mkdir -p $outdir/$sampid/0.shell/a05.mutect2") unless (-e "$outdir/$sampid/0.shell/a05.mutect2");
        system("mkdir -p $outdir/$sampid/0.shell/a06.strelka")unless(-e "$outdir/$sampid/0.shell/a06.strelka");
        system("mkdir -p $outdir/$sampid/0.shell/a07.strelka2") unless (-e "$outdir/$sampid/0.shell/a07.strelka2");
        system("mkdir -p $outdir/$sampid/0.shell/a08.svaba") unless (-e "$outdir/$sampid/0.shell/a08.svaba");
        system("mkdir -p $outdir/$sampid/0.shell/a09.muse") unless (-e "$outdir/$sampid/0.shell/a09.muse");
        system("mkdir -p $outdir/$sampid/0.shell/a10.mergeMut") unless (-e "$outdir/$sampid/0.shell/a10.mergeMut");
        system("mkdir -p $outdir/$sampid/0.shell/a11.hlasomatic") unless (-e "$outdir/$sampid/0.shell/a11.hlasomatic");
        system("mkdir -p $outdir/$sampid/0.shell/a12.kraken2") unless (-e "$outdir/$sampid/0.shell/a12.kraken2");
        ## somatic mutation
        system("mkdir -p $outdir/$sampid/3.somatic/mutect") unless (-e "$outdir/$sampid/3.somatic/mutect");
        system("mkdir -p $outdir/$sampid/3.somatic/mutect2") unless (-e "$outdir/$sampid/3.somatic/mutect2");
        system("mkdir -p $outdir/$sampid/3.somatic/strelka2") unless (-e "$outdir/$sampid/3.somatic/strelka2");
        system("mkdir -p $outdir/$sampid/3.somatic/muse") unless (-e "$outdir/$sampid/3.somatic/muse");
        system("mkdir -p $outdir/$sampid/3.somatic/result") unless (-e "$outdir/$sampid/3.somatic/result");
        system("mkdir -p $outdir/$sampid/3.somatic/hlasomatic") unless (-e "$outdir/$sampid/3.somatic/hlasomatic");   
        ##qc
        system("mkdir -p $outdir/$sampid/2.alignment/QC") unless (-e "$outdir/$sampid/2.alignment/QC");
    }
}

sub generateShell{
    my ($output_shell, $content, $finish_string) = @_;
    my $shell_name = basename($output_shell);
    $finish_string ||= "Still_waters_run_deep";
    open OUT, ">$output_shell" or die "Cannot open file $output_shell:$!";
    print OUT "#!/bin/bash\n";
    print OUT "echo hostname: `hostname`\n";
    print OUT "echo ==========start at : `date` ==========\n";
    print OUT "$content && \\\n";
    print OUT "echo ==========end at : `date` ========== && \\\n";
    print OUT "echo $finish_string 1>&2 && \\\n";
    print OUT "echo $finish_string > $output_shell.sign\n";
    # print OUT "qstat -j $shell_name | grep \"usage\" > $output_shell.log";
    print OUT "qstat -j $shell_name > $output_shell.log";    
    close OUT;
}

sub runsplit{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        for my $dtype(sort keys %{$hash{$sampid}}){
            for my $lane(sort keys %{$hash{$sampid}{$dtype}}){
                my @fq1 = @{$hash{$sampid}{$dtype}{$lane}{fq1}};
                my @fq2 = @{$hash{$sampid}{$dtype}{$lane}{fq2}};
                for my $i(0..$#fq1){
                    my $tmp=substr($fq1[$i],0,-3);my $readsnum;
                    if(-e "$tmp.fqStat.txt"){
                        $readsnum=`cat $tmp.fqStat.txt | grep 'ReadNum' | cut -f2`;chomp($readsnum);
                        $splitfqs{$sampid}{$dtype}{$lane}{fqstat_exist} = "true";
                    }else{
                        push @{$splitfqs{$sampid}{$dtype}{$lane}{fq1}}, $fq1[$i];
                        push @{$splitfqs{$sampid}{$dtype}{$lane}{fq2}}, $fq2[$i];
                        $splitfqs{$sampid}{$dtype}{$lane}{fqstat_exist} = "false";
                        next;
                    }
                    my $splitrow=POSIX::ceil ($readsnum/$configs{split_num});
                    $splitrow *=4;
                    my ($unzipfile, $unzipfile1);
                    if(($splitrow*($configs{split_num}-1)/4)<$readsnum){
                        my $shell="$outdir/$sampid/0.shell/a01.clean/a01.fqsplit_$sampid\_$dtype\_$lane\_$i.1.sh";
                        my $cmd;
                        if($fq1[$i]=~/\.gz/){
                            $unzipfile=basename $fq1[$i];
                            $unzipfile=substr($unzipfile,0,-3);
                            $cmd = "gunzip -c $fq1[$i] > $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile && \\\n";
                            $cmd.= "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile -d -a 2";
                            $cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_ && \\\n";
                            $cmd.= "rm -rf $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile";
                        }else{
                            $unzipfile=$fq1[$i];
                            $cmd = "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile -d -a 2";
                            $cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_";
                        }
                        generateShell($shell,$cmd);
                        $scripts{split}{$sampid}{$dtype}{$lane}{$i.1}=$shell;
                        $shell="$outdir/$sampid/0.shell/a01.clean/a01.fqsplit_$sampid\_$dtype\_$lane\_$i.2.sh";
                        undef $cmd;
                        if($fq2[$i]){
                            if($fq2[$i]=~/\.gz/){
                                $unzipfile1=basename $fq2[$i];
                                $unzipfile1=substr($unzipfile1,0,-3);
                                $cmd = "gunzip -c $fq2[$i] > $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1 && \\\n";
                                $cmd.= "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1 -d -a 2";
                                $cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_ && \\\n";
                                $cmd.= "rm -rf $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1";
                            }else{
                                $unzipfile1=$fq2[$i];
                                $cmd = "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1 -d -a 2";
                                $cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_";
                            }
                        }
                        generateShell($shell,$cmd);
                        $scripts{split}{$sampid}{$dtype}{$lane}{$i.2}=$shell;
                    }else{
                        print "split err: the split num is too high for $fq1[$i]";
                        print " and $fq2[$i]" if($fq2[$1]);
                        print "\n";
                    }
                    my $s = $configs{split_num} -1;
                    for my $n(0..$s){
                        my ($splitfq1,$splitfq2);
                        if($n < 10){
                            $splitfq1 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_0$n";
                            $splitfq2 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_0$n" if($unzipfile1);
                        }else{
                            $splitfq1 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_$n";
                            $splitfq2 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_$n" if($unzipfile1);
                        }
                        push @{$splitfqs{$sampid}{$dtype}{$lane}{fq1}}, $splitfq1;
                        push @{$splitfqs{$sampid}{$dtype}{$lane}{fq2}}, $splitfq2;
                    }
                }
            }
        }
    }
}

sub FastqClean{
    my (%hash) = @_;
    # my %hash = %{$hash};
    for my $sampid(sort keys %hash){
        for my $dtype(sort keys %{$hash{$sampid}}){
            # $cleanfqs{$sampid}{$dtype} = "$outdir/$sampid/2.alignment/$sampid.$dtype.manifest.csv";
            # open ELIST, ">$outdir/$sampid/2.alignment/$sampid.$dtype.manifest.csv" or die "$!";
            for my $lane(sort keys %{$hash{$sampid}{$dtype}}){
                my @fq1 = @{$hash{$sampid}{$dtype}{$lane}{fq1}};
                my @fq2 = @{$hash{$sampid}{$dtype}{$lane}{fq2}};
                for my $i(0..$#fq1){
                    my ($out_fq1, $out_fq2);
                    my $mydir = "$outdir/$sampid/1.clean/$dtype/$lane/fq$i";
                    system("mkdir -p $mydir") unless (-e "$mydir");
                    $out_fq1 = "$mydir/clean_1.$i.fq.gz";
                    $out_fq2 = "$mydir/clean_2.$i.fq.gz";
                    my $shell = "$outdir/$sampid/0.shell/a01.clean/a01.cleanfq.$sampid\_$dtype\_$lane\_$i.sh";
                    my $cmd;
                    if($configs{soapnukeclean} =~ /false/){
                        $cmd.= "$configs{fastp} -i $fq1[$i] -o $out_fq1";
                        $cmd.= " -I $fq2[$i] -O $out_fq2" if($fq2[$i]);
                        $cmd.= " -j $mydir/$sampid.$dtype.$lane.fastp.json";
                        $cmd.= " -h $mydir/$sampid.$dtype.$lane.fastp.html";
                        $cmd.= " $configs{fastpparameter}";
                    }else{
                        $cmd.= "$configs{soapnuketool} filter -1 $fq1[$i] -C clean_1.$i.fq.gz";
                        $cmd.= " -2 $fq2[$i] -D clean_2.$i.fq.gz" if($fq2[$i]);
                        $cmd.= " $configs{soapnukeparameter} -o $mydir";
                    }
                    if($configs{split_num} > 1 && $hash{$sampid}{$dtype}{$lane}{fqstat_exist} =~ /true/i){
                        $cmd.= "&& \\\nrm -rf $fq1[$i]";
                        $cmd.= " && \\\nrm -rf $fq2[$i]" if($fq2[$i]);

                    }
                    generateShell($shell,$cmd);
                    $scripts{clean}{$sampid}{$dtype}{$lane}{$i}=$shell;
                    push @{$cleanfastq{$sampid}{$dtype}{fq1}}, $out_fq1;
                    push @{$cleanfastq{$sampid}{$dtype}{fq2}}, $out_fq2;
                    # print ELIST 
                }
            }
            # my $mega_fq1 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq1}});
            # my $mega_fq2 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq2}});
            # print ELIST join("\t","$sampid\-$dtype",$mega_fq1);
            # print ELIST "\t$mega_fq2" unless($mega_fq2 =~ /,{2,}/);
            # print ELIST "\n";
            # close ELIST;
        }
    }
}

sub FastqLink{
    my ($hash) = @_;
    my %hash = %{$hash};
    for my $sampid(sort keys %hash){
        for my $dtype(sort keys %{$hash{$sampid}}){
            # if ($dtype =~ /Tumor/i or $dtype =~ /Normal/i){
            #     $cleanfqs{$sampid}{$dtype} = "$outdir/$sampid/2.alignment/$sampid.$dtype.manifest.csv";
            #     open ELIST, ">$outdir/$sampid/2.alignment/$sampid.$dtype.manifest.csv" or die "$!";
            # }
            for my $lane(sort keys %{$hash{$sampid}{$dtype}}){
                my @fq1 = @{$hash{$sampid}{$dtype}{$lane}{fq1}};
                my @fq2 = @{$hash{$sampid}{$dtype}{$lane}{fq2}};
                for my $i(0..$#fq1){
                    my ($out_fq1, $out_fq2);
                    $out_fq1 = "$outdir/$sampid/1.clean/$dtype/$lane/clean_1.$i.fq.gz";
                    $out_fq2 = "$outdir/$sampid/1.clean/$dtype/$lane/clean_2.$i.fq.gz" if($fq2[$i]);
                    `ln -s $fq1[$i] $out_fq1`;
                    `ln -s $fq2[$i] $out_fq2` if($fq2[$i]);
                    push @{$cleanfastq{$sampid}{$dtype}{fq1}}, $out_fq1;
                    push @{$cleanfastq{$sampid}{$dtype}{fq2}}, $out_fq2;
                }
            }
            # my $mega_fq1 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq1}});
            # my $mega_fq2 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq2}});
            # print ELIST join("\t","$sampid\-$dtype",$mega_fq1);
            # print ELIST "\t$mega_fq2" unless($mega_fq2 =~ /,{2,}/);
            # print ELIST "\n";
            # close ELIST;
        }
    }
}

sub FastqCorrect{
    # only suit for PE reads
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        for my $dtype(keys %{$hash{$sampid}}){
            my (@fq1, @fq2);
            my $mydir = "$outdir/$sampid/1.clean/$dtype/Cor";
            system("mkdir -p $mydir") unless(-e $mydir);
            @fq1 = @{$hash{$sampid}{$dtype}{fq1}};
            @fq2 = @{$hash{$sampid}{$dtype}{fq2}};
            for my $i(0..$#fq1){
                my $shell = "$outdir/$sampid/0.shell/a01.clean/a01.cleanfq.$sampid\_$dtype\_$i.sh";
                my $cmd = "if [ -e $fq2[$i] ]\nthen\n";
                $cmd.= "\tperl $configs{rcorrector} -t 8 -1 $fq1[$i] -2 $fq2[$i] -od $mydir && \\\n";
                $cmd.= "\t$configs{python} $configs{filterfq} -1 $mydir/clean_1.$i.cor.fq.gz -2 $mydir/clean_2.$i.cor.fq.gz -s ${sampid}_${dtype}_$i && \\\n";
                $cmd.= "\tgzip $mydir/unfixrm_clean_1.$i.cor.fq && \\\n";
                $cmd.= "\tgzip $mydir/unfixrm_clean_2.$i.cor.fq && \\\n";
                $cmd.= "\trm $fq1[$i] $fq2[$i]\n";
                $cmd.= "fi";
                my $out_fq1 = "$mydir/unfixrm_clean_1.$i.cor.fq.gz";
                my $out_fq2 = "$mydir/unfixrm_clean_2.$i.cor.fq.gz";
                push @{$correctfastq{$sampid}{$dtype}{fq1}}, $out_fq1;
                push @{$correctfastq{$sampid}{$dtype}{fq2}}, $out_fq2;
                generateShell($shell, $cmd);
                $scripts{correct}{$sampid}{$dtype}{$i}=$shell;
            }
        }
    }
}

sub Alignment{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        for my $dtype(keys %{$hash{$sampid}}){
            my $mydir = "$outdir/$sampid/2.alignment/$dtype";
            system("mkdir -p $mydir/tmp") unless(-e "$mydir/tmp");
            my $shell = "$outdir/$sampid/0.shell/a02.alignment/a02.alignment.$sampid\_$dtype.sh";
            my $cmd;
            my $bam;
            if($configs{alignment} =~ /true/i){
                # run STAR alignment for fastq input
                my (@clean_fq1, @clean_fq2);
                @clean_fq1 = @{$hash{$sampid}{$dtype}{fq1}};
                @clean_fq2 = @{$hash{$sampid}{$dtype}{fq2}};
                $cleanfqs{$sampid}{$dtype} = "$outdir/$sampid/2.alignment/$sampid.$dtype.manifest.csv";
                open ELIST, "> $cleanfqs{$sampid}{$dtype}" or die "$!";
                for my $i(0..$#clean_fq1){
                    my ($rgid, $rgsm, $rgpl, $rgpu, $rglb) = ("$sampid\_$dtype\_fq$i", "$sampid\_$dtype", "ILLUMINA", "unit$i", "unit$i");
                    my $fq1 = $clean_fq1[$i];
                    my $fq2 = $clean_fq2[$i] || "-";
                    print ELIST "$fq1\t$fq2\tID:$rgid\tSM:$rgsm\tPL:$rgpl\tPU:$rgpu\tLB:$rglb\n";
                }
                close ELIST;
                # $cmd.= "$configs{star} --runThreadN 12 --runMode alignReads --genomeDir $configs{star_ref} --outFileNamePrefix $mydir/$sampid.$dtype. --outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 100 --readFilesCommand \"gunzip -c\" $configs{star_par} --readFilesManifest $cleanfqs{$sampid}{$dtype} && \\\n";
                $cmd.= "$configs{star} --runThreadN 12 --runMode alignReads --genomeDir $configs{star_ref} --outFileNamePrefix $mydir/$sampid.$dtype. --outSAMtype BAM Unsorted --readFilesCommand \"gunzip -c\" $configs{star_par} --readFilesManifest $cleanfqs{$sampid}{$dtype} --outSAMattributes NH HI AS nM RG && \\\n";
                $cmd.= "$configs{samtools} sort -\@ 8 -o $mydir/$sampid.$dtype.Aligned.sortedByCoord.out.bam $mydir/$sampid.$dtype.Aligned.out.bam && \\\n";
                $cmd.= "$configs{samtools} index -b -\@ 8 $mydir/$sampid.$dtype.Aligned.sortedByCoord.out.bam && \\\n";
                $cmd.= "gzip $mydir/$sampid.$dtype.Unmapped.out.mate* && \\\n";
                $cmd.= "rm $mydir/$sampid.$dtype.Aligned.out.bam";
                $bam = "$mydir/$sampid.$dtype.Aligned.sortedByCoord.out.bam";
            }else{
                # for BAM input
                my @this_bams = @{$hash{$sampid}{$dtype}{bam}};
                # rewrite RG contents
                for my $i(0..$#this_bams){
                    my ($rgid, $rgsm, $rgpl, $rgpu, $rglb, $this_bam) = ("$sampid\_$dtype\_bam$i", "$sampid\_$dtype", "ILLUMINA", "unit$i", "unit$i", "$this_bams[$i]");
                    if($configs{addreadgroup} =~ /true/){
                        $cmd.= "$configs{java} -jar $configs{picard} AddOrReplaceReadGroups I=$this_bam O=$mydir/tmp/sort.$i.bam RGID=$rgid RGSM=$rgsm RGPL=$rgpl RGPU=$rgpu RGLB=$rglb SORT_ORDER=coordinate && \\\n";
                    }else{
                        $cmd.= "ln -s $this_bam $mydir/tmp/sort.$i.bam && \\\n"
                    }
                    $cmd.= "$configs{samtools} index -b -\@ 8 $mydir/tmp/sort.$i.bam && \\\n";
                }
                # all BAM files merging, sorting and indexing
                $cmd.= "$configs{samtools} merge -f -\@ 8 --write-index $mydir/$sampid.$dtype.Aligned.sortedByCoord.out.bam";
                for my $i(0..$#this_bams){
                    $cmd .= " $this_bams[$i]";
                }
                $cmd.= " && \\\n";
                # $cmd.= "$configs{samtools} sort -\@ 8 -T $mydir/tmp/sampid.$dtype -o $mydir/$sampid.$dtype.sortedByCoord.out.bam $mydir/$sampid.$dtype.merged.bam && \\\n";
                # $cmd.= "$configs{samtools} index -\@ 8 $mydir/$sampid.$dtype.sortedByCoord.out.bam";
                $cmd.= "rm $mydir/$sampid.$dtype.merged.bam  $mydir/tmp/*";
                $bam = "$mydir/$sampid.$dtype.Aligned.sortedByCoord.out.bam";
            }
            $bqsrbam{$sampid}{$dtype} = $bam;
            $unmapfqs{$sampid}{$dtype}{fq1} = "$mydir/$sampid.$dtype.Unmapped.out.mate1.gz";
            $unmapfqs{$sampid}{$dtype}{fq2} = "$mydir/$sampid.$dtype.Unmapped.out.mate2.gz";
            $scripts{alignment}{$sampid}{$dtype} = $shell;
            generateShell($shell,$cmd);
            if($configs{bampostprocess} =~ /true/i){
                for my $mychr(@chrs){
                    $shell = "$outdir/$sampid/0.shell/a02.alignment/a02.alignment.postprocessing.$sampid\_$dtype.$mychr.sh";
                    $cmd = "$configs{samtools} view -\@ 4 -bh -o $mydir/$sampid.$dtype.$mychr.sort.bam -L $configs{reference_block_bed_path}/$mychr.bed --write-index $mydir/$sampid.$dtype.Aligned.sortedByCoord.out.bam && \\\n";
                    # Mark Duplicates
                    $cmd.= "$configs{GATK4} MarkDuplicates --INPUT $mydir/$sampid.$dtype.$mychr.sort.bam --OUTPUT $mydir/$sampid.$dtype.$mychr.sort.dedup.bam --CREATE_INDEX true --TMP_DIR $mydir/tmp --VALIDATION_STRINGENCY SILENT --METRICS_FILE $mydir/$sampid.$dtype.$mychr.dedupped.metrics && \\\n";
                    $cmd.= "rm $mydir/$sampid.$dtype.$mychr.sort.bam*\n";
                    # Split N Cigar Reads
                    $cmd.= "$configs{GATK4} --java-options \"-Xms10G -Djava.io.tmpdir=$mydir/tmp -XX:-UseGCOverheadLimit -XX:GCHeapFreeLimit=10\" SplitNCigarReads -R $configs{reference} -I $mydir/$sampid.$dtype.$mychr.sort.dedup.bam -O $mydir/$sampid.$dtype.$mychr.sort.dedup.splitN.bam -L $configs{reference_block_bed_path}/$mychr.bed --tmp-dir $mydir/tmp && \\\n";
                    $cmd.= "rm $mydir/$sampid.$dtype.$mychr.sort.dedup.bam*\n";
                    # BaseRecalibrator
                    $cmd.= "$configs{GATK4} --java-options \"-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -Xms4G\" BaseRecalibrator -R $configs{reference} -I $mydir/$sampid.$dtype.$mychr.sort.dedup.splitN.bam -L $configs{reference_block_bed_path}/$mychr.bed -O $mydir/$sampid.$dtype.$mychr.recal_data.csv";
                    $cmd.= " --known-sites $configs{Cosmic}" if($configs{Cosmic});
                    $cmd.= " --known-sites $configs{GATKdbsnp}" if($configs{GATKdbsnp});
                    $cmd.= " --known-sites $configs{'1kg_phase1_snp'}" if($configs{'1kg_phase1_snp'});
                    $cmd.= " --known-sites $configs{Mills_indel}" if($configs{Mills_indel});
                    $cmd.= " && \\\n";
                    # ApplyBQSR
                    $cmd.= "$configs{GATK4} --java-options \"-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3G\" ApplyBQSR --add-output-sam-program-record -R $configs{reference} -I $mydir/$sampid.$dtype.$mychr.sort.dedup.splitN.bam -L $configs{reference_block_bed_path}/$mychr.bed -O $mydir/$sampid.$dtype.$mychr.sort.dedup.splitN.bqsr.bam --bqsr-recal-file $mydir/$sampid.$dtype.$mychr.recal_data.csv && \\\n";
                    $cmd.= "rm $mydir/$sampid.$dtype.$mychr.sort.dedup.splitN.bam*";
                    $scripts{bampostprocess}{$sampid}{$dtype}{$mychr} = $shell;
                    $SplitedBams{$sampid}{$dtype}{$mychr} = "$mydir/$sampid.$dtype.$mychr.sort.dedup.splitN.bqsr.bam";
                    generateShell($shell, $cmd);
                }
                # merge bqst BAM files
                $shell = "$outdir/$sampid/0.shell/a02.alignment/a02.bqsrbamMerge.$sampid\_$dtype.sh";
                $cmd = "$configs{samtools} merge -f -\@ 8 --write-index $mydir/$sampid.$dtype.RNA.sort.dedup.splitN.bqsr.final.bam";
                for my $mychr(@chrs){
                    $cmd.= " $SplitedBams{$sampid}{$dtype}{$mychr}";
                }
                $cmd.= " && \\\n";
                $cmd.= "$configs{samtools} index -b -\@ 8 $mydir/$sampid.$dtype.RNA.sort.dedup.splitN.bqsr.final.bam && \\\n";
                $cmd.= "rm $mydir/$sampid.$dtype.*.sort.dedup.splitN.bqsr.bam*";
                $bqsrbam{$sampid}{$dtype} = "$mydir/$sampid.$dtype.RNA.sort.dedup.splitN.bqsr.final.bam";
                $scripts{BamMerge}{$sampid}{$dtype} = $shell;
                generateShell($shell, $cmd);
                # merge bqsr reports
                $shell = "$outdir/$sampid/0.shell/a02.alignment/a02.BqsrReportMerge.$sampid\_$dtype.sh";
                $cmd = "$configs{GATK4} --java-options \"-Xms4G\" GatherBQSRReports";
                for my $mychr(@chrs){
                    $cmd.= " -I $mydir/$sampid.$dtype.$mychr.recal_data.csv";
                }
                $cmd.= " -O $mydir/$sampid.$dtype.recal_data.csv && \\\n";
                $cmd.= "rm $mydir/$sampid.$dtype.*.recal_data.csv";
                $scripts{GatherBQSR}{$sampid}{$dtype} = $shell;
                generateShell($shell, $cmd);
            }
        }
    }
}

sub bamqc{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/2.alignment/QC";
        if($configs{qcvcf}){
            my $tbam = $hash{$sampid}{Tumor};
            my $nbam = $hash{$sampid}{Normal};
            system("mkdir -p $mydir/tmp") unless (-e "$mydir/tmp");
            my $shell = "$outdir/$sampid/0.shell/a02.alignment/QC.sh";
            my $cmd = "$configs{python} $configs{qctool} --bam1 $nbam --bam2 $tbam --output $mydir/bam_matcher.report.txt --vcf $configs{qcvcf} --reference $configs{reference} --scratch-dir $mydir/tmp $configs{qcpar} && \\\n";
            $cmd.= "rm -rf $mydir/tmp";
            generateShell($shell,$cmd);
            $scripts{BAMqc}{$sampid} = $shell;
        }
        for my $dtype(sort keys %{$hash{$sampid}}){
            my $shell = "$outdir/$sampid/0.shell/a02.alignment/$dtype.coverage.sh";
            my $cmd = "$configs{samtools} bedcov $configs{targetRegion} $hash{$sampid}{$dtype} |gzip -c > $mydir/$sampid.$dtype.bedcov.txt.gz && \\\n";
            $cmd.= "$configs{samtools} stats -\@ 4 $hash{$sampid}{$dtype} > $mydir/$sampid.$dtype.coverage.txt.gz";
            generateShell($shell,$cmd);
            $scripts{coverage}{$sampid}{$dtype} = $shell            
        }
    }
}

sub svaba{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/svaba";
        system("mkdir -p $mydir") unless (-e "$mydir");
        my $tbam = $hash{$sampid}{Tumor};
        my $nbam = $hash{$sampid}{Normal};
        for my $mychr(@chrs){
            my $shell = "$outdir/$sampid/0.shell/a08.svaba/svaba.$mychr.sh";
            my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n"; # [Ming-20191129] source
            $cmd.= "$configs{svabatool} run -G $configs{reference} -t $tbam -n $nbam -k $configs{reference_block_bed_path}/$mychr.bed -a $mydir/$mychr";
            $cmd.= " -D $configs{GATKdbsnp}" if($configs{GATKdbsnp});
            generateShell($shell,$cmd);
            $scripts{svaba}{$sampid}{$mychr} = $shell;
        }
        my $shell = "$outdir/$sampid/0.shell/a08.svaba/merge.svaba.sh";
        my $cmd = "less $mydir/chr1.svaba.somatic.indel.vcf | grep '^##' | grep -v '^##source'> $mydir/svaba.head && \\\n";
        $cmd.= "less $mydir/chr*.svaba.somatic.indel.vcf | grep '^##source' >> $mydir/svaba.head && \\\n";
        $cmd.= "less $mydir/chr1.svaba.somatic.indel.vcf | grep '^#CHROM' >> $mydir/svaba.head && \\\n";
        $cmd.= "less $mydir/chr*.svaba.somatic.indel.vcf | grep -v '#' | grep PASS | cat $mydir/svaba.head - > $mydir/$sampid.svaba.vcf && \\\n";
        $cmd.= "rm $mydir/chr*";
        generateShell($shell,$cmd);
        $scripts{svaba}{$sampid}{merge} = $shell;
        $indelvcf{$sampid}{svaba}="$mydir/$sampid.svaba.vcf";
    }
}

sub mutectSNV{
    my (%hash) = @_;    # BQSR bam, whole Tumor/Normal bam
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/mutect";
        my $tbam = $hash{$sampid}{Tumor};
        my $nbam = $hash{$sampid}{Normal};
        for my $mychr(@chrs){
            my $shell = "$outdir/$sampid/0.shell/a04.mutect/mutect.bqsr.$mychr.sh";
            my $thistmp = "$mydir/$mychr.tmp";
            system("mkdir -p $thistmp") unless (-e "$thistmp");
            my $cmd = "rm -f $mydir/$sampid.mutect.$mychr.txt.gz\n";
            $cmd.= "$configs{java17} -Xmx10g -Djava.io.tmpdir=$thistmp -XX:-UseGCOverheadLimit -jar $configs{mutectjar} -T MuTect -R $configs{reference} --input_file:normal $nbam --input_file:tumor $tbam --normal_sample_name $sampid\_Normal --tumor_sample_name $sampid\_Tumor --vcf $mydir/$sampid.mutect.$mychr.vcf.gz --out $mydir/$sampid.mutect.$mychr.txt";
            $cmd.= " --dbsnp $configs{dbsnp}" if($configs{dbsnp});
            $cmd.= " --cosmic $configs{Cosmic}" if($configs{Cosmic});
            $cmd.= " -L $configs{reference_block_bed_path}/$mychr.bed --enable_extended_output --downsampling_type NONE --filter_reads_with_N_cigar && \\\n";
            $cmd.= "gzip $mydir/$sampid.mutect.$mychr.txt";
            generateShell($shell,$cmd);
            $scripts{mutect}{step1}{$sampid}{$mychr} = $shell;
        }
        my $shell = "$outdir/$sampid/0.shell/a04.mutect/mergevcf.mutect.sh";
        my $cmd = "less $mydir/$sampid.mutect.chr1.txt.gz | head -2 > $mydir/head.txt && \\\n";
        $cmd.= "less $mydir/$sampid.mutect.*.txt.gz | grep \"$sampid\" | cat $mydir/head.txt - | gzip -f > $mydir/MuTect.txt.gz && \\\n";
        $cmd.= "less $mydir/$sampid.mutect.chr1.vcf.gz | grep \"#\" > $mydir/head2.txt && \\\n";
        $cmd.= "less $mydir/$sampid.mutect.*.vcf.gz | grep -v \"#\" | cat $mydir/head2.txt - | gzip -f > $mydir/$sampid.MuTect.vcf.gz && \\\n";
        $cmd.= "rm -rf $mydir/*.tmp && \\\n";
        $cmd.= "rm $mydir/head* && \\\n";
        $cmd.= "rm $mydir/$sampid.mutect.*.txt.gz && \\\n";
        $cmd.= "rm $mydir/$sampid.mutect.*.vcf.gz && \\\n";
        $cmd.= "rm $mydir/$sampid.mutect.*.vcf.gz.idx";
        generateShell($shell,$cmd);
        $scripts{mutect}{step2}{$sampid} = $shell;
        $snvvcf{$sampid}{mutect} = "$mydir/$sampid.MuTect.vcf.gz";
    }
}

sub mutect2{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/mutect2";
        system("mkdir -p $mydir/contamination/tmp") unless (-e "$mydir/contamination/tmp");
        system("mkdir -p $mydir/artifact/tmp") unless (-e "$mydir/artifact/tmp");
        system("mkdir -p $mydir/PoN/tmp") unless (-e "$mydir/PoN/tmp");
        system("mkdir -p $mydir/tmp") unless (-e "$mydir/tmp");
        my $tbam = $hash{$sampid}{Tumor};
        my $nbam = $hash{$sampid}{Normal};
        ## Contamination Estimate
        my $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.contamination.sh";
        my $cmd = "#$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\" GetPileupSummaries -I $tbam -V $configs{GetPileupSummaries} -L $configs{reference_bed} -O $mydir/contamination/$sampid\_Tumor.pileups.table && \\\n";
        $cmd.= "#$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\" GetPileupSummaries -I $nbam -V $configs{GetPileupSummaries} -L $configs{reference_bed} -O $mydir/contamination/$sampid\_Normal.pileups.table && \\\n";
        $cmd.= "#$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\" CalculateContamination -I $mydir/contamination/$sampid\_Tumor.pileups.table -matched $mydir/contamination/$sampid\_Normal.pileups.table -O $mydir/contamination/$sampid.contamination.table";
        generateShell($shell,$cmd);
        $scripts{mutect2}{$sampid}{contamination} = $shell;
        ## Artifact calculation
        # $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.artifact.sh";
        # $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        # $cmd.= "$configs{GATK4} --java-options \"-Xmx4g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/artifact/tmp\" CollectSequencingArtifactMetrics -R $configs{reference} -I $tbam --FILE_EXTENSION \".txt\" -O $mydir/artifact/$sampid.artifact";
        # generateShell($shell,$cmd);
        # $scripts{mutect2}{$sampid}{artifact} = $shell;
        ##normal Panel
        # $shell = "$outdir/$sampid/0.shell/a06.mutect2/CreateSomaticPanelOfNormals.sh";
        # $cmd = "mkdir -p $mydir/PoN/tmp\n";
        # $cmd.= "$configs{GATK4} --java-options \"-Xmx16g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5";
        # $cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/PoN/tmp\" Mutect2 -R $configs{reference}";
        # $cmd.= " -I $nbam -tumor $sampid\_Normal";
        # $cmd.= " --germline-resource $configs{gnomAD}" if($configs{gnomAD});
        # $cmd.= " -O $mydir/PoN/$sampid\_Normal.vcf.gz -L $configs{targetRegion}";
        # generateShell($shell,$cmd) unless($configs{panel_of_normal});
        # $scripts{mutect2}{$sampid}{PoN} = $shell;
        ## Mutect2
        my @vcfs;
        for my $mychr(@chrs){
            my $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.call.$mychr.sh";
            my $thistmp = "$mydir/$mychr.tmp";
            system("mkdir -p $thistmp") unless (-e "$thistmp");
            $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
            # CalculateContamination
            $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$thistmp\" GetPileupSummaries -I $tbam -V $configs{GetPileupSummaries} -L $configs{reference_block_bed_path}/$mychr.bed -O $mydir/contamination/$sampid\_Tumor.$mychr.pileups.table && \\\n";
            $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$thistmp\" GetPileupSummaries -I $nbam -V $configs{GetPileupSummaries} -L $configs{reference_block_bed_path}/$mychr.bed -O $mydir/contamination/$sampid\_Normal.$mychr.pileups.table && \\\n";
            $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$thistmp\" CalculateContamination -I $mydir/contamination/$sampid\_Tumor.$mychr.pileups.table -matched $mydir/contamination/$sampid\_Normal.$mychr.pileups.table -O $mydir/contamination/$sampid.contamination.$mychr.table && \\\n";
            # Mutect2 call mutation
            $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$thistmp\" Mutect2 -R $configs{reference} -I $tbam -tumor $sampid\_Tumor -I $nbam -normal $sampid\_Normal";
            $cmd.= " --germline-resource $configs{gnomAD}" if($configs{gnomAD});
            if($configs{panel_of_normal}){
                $cmd.= " --panel-of-normals $configs{panel_of_normal}";
            }
            $cmd.= " -L $configs{reference_block_bed_path}/$mychr.bed -O $mydir/$sampid.mutect2.$mychr.vcf.gz && \\\n";
            # FilterMutectCalls
            $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$thistmp\" FilterMutectCalls -V $mydir/$sampid.mutect2.$mychr.vcf.gz -R $configs{reference} -contamination-table $mydir/contamination/$sampid.contamination.$mychr.table -O $mydir/$sampid.mutect2.$mychr.filter1.vcf.gz";
            generateShell($shell,$cmd);
            $scripts{mutect2}{step1}{$sampid}{$mychr} = $shell;
            push @vcfs, "$mydir/$sampid.mutect2.$mychr.filter1.vcf.gz";
        }
        $shell = "$outdir/$sampid/0.shell/a05.mutect2/mergevcf.mutect2.sh";
        $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n"; # [Ming-20191203] source
        $cmd.= "$configs{vcfconcat} @vcfs | $configs{vcfsort} -t $mydir > $mydir/$sampid.mutect2.filter1.vcf && \\\n";
        $cmd.= "$configs{bgzip} $mydir/$sampid.mutect2.filter1.vcf && \\\n";
        $cmd.= "$configs{tabix} $mydir/$sampid.mutect2.filter1.vcf.gz";
        generateShell($shell,$cmd);
        $scripts{mutect2}{step2}{$sampid} = $shell;
        ## Filter Mutect Calls
        $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.filter.sh";
        # $cmd = "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterMutectCalls -V $mydir/$sampid.mutect2.vcf.gz -contamination-table $mydir/contamination/$sampid.contamination.table -O $mydir/$sampid.mutect2.filter1.vcf.gz && \\\n";
        $cmd = "#$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterMutectCalls -V $mydir/$sampid.mutect2.vcf.gz -contamination-table $mydir/contamination/$sampid.contamination.table -O $mydir/$sampid.mutect2.filter2.vcf.gz && \\\n";
        # $cmd.= "$configs{GATK4} --java-options \"-Xmx50g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterByOrientationBias --artifact-modes \"G/T\" --artifact-modes \"C/T\" -V $mydir/$sampid.mutect2.filter1.vcf.gz -P $mydir/artifact/$sampid.artifact.pre_adapter_detail_metrics.txt -O $mydir/$sampid.mutect2.filter2.vcf.gz && \\\n";
        $cmd.= "rm -rf $mydir/chr*.tmp && \\\nrm $mydir/$sampid.mutect2.chr* ";
        generateShell($shell,$cmd);
        $scripts{mutect2}{$sampid}{filter} = $shell;
        $snvvcf{$sampid}{mutect2} = "$mydir/$sampid.mutect2.filter1.vcf.gz";
        $indelvcf{$sampid}{mutect2} = "$mydir/$sampid.mutect2.filter1.vcf.gz";
    }
}

sub strelka2{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/strelka2";
        my $tbam = $hash{$sampid}{Tumor};
        my $nbam = $hash{$sampid}{Normal};
        my $shell = "$outdir/$sampid/0.shell/a07.strelka2/$sampid.strelka2.sh";
        my $cmd = "rm -rf $mydir/*\n";
        $cmd.= "$configs{python} $configs{strelka2py} --normalBam $nbam --tumorBam $tbam --referenceFasta $configs{reference} --runDir $mydir && \\\n";
        $cmd.= "$configs{python} $mydir/runWorkflow.py $configs{runstrelka2flow}";
        generateShell($shell,$cmd);
        $scripts{strelka2}{$sampid} = $shell;
        $indelvcf{$sampid}{strelka2} = "$mydir/results/variants/somatic.indels.vcf.gz";
        $snvvcf{$sampid}{strelka2} = "$mydir/results/variants/somatic.snvs.vcf.gz";
    }
}

sub muse{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/muse";
        my $tbam = $hash{$sampid}{Tumor};
        my $nbam = $hash{$sampid}{Normal};
        for my $mychr(@chrs){
            my $shell="$outdir/$sampid/0.shell/a09.muse/muse.$mychr.sh";
            my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n"; # [Ming-20191129] source
            $cmd.= "$configs{muse} call -f $configs{reference} -l $configs{reference_block_bed_path}/$mychr.bed -O $mydir/$mychr $tbam $nbam";
            # $cmd.= "$configs{muse} sump -I $mydir/$mychr.MuSE.txt -G -O $mydir/$mychr.vcf -D $configs{GATKdbsnp}";
            generateShell($shell,$cmd);
            $scripts{MuSE}{step1}{$sampid}{$mychr}=$shell;
        }
        my $shell = "$outdir/$sampid/0.shell/a09.muse/mergemuse.sh";
        my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        $cmd.= "less $mydir/chr1.MuSE.txt | grep '^#'> $mydir/txt.head.txt && \\\n";
        $cmd.= "cat $mydir/chr*.MuSE.txt | grep -v '^#' | cat $mydir/txt.head.txt - > $mydir/all.Muse.txt && \\\n";
        $cmd.= "$configs{muse} sump -I $mydir/all.Muse.txt -G -O $mydir/all.MuSE.vcf -D $configs{GATKdbsnp} && \\\n";
        $cmd.= "less $mydir/all.MuSE.vcf | grep '^#' > $mydir/head.txt && \\\n";
        $cmd.= "cat $mydir/all.MuSE.vcf | grep -v '^#' | grep 'PASS'| cat $mydir/head.txt - | gzip -f > $mydir/$sampid.MuSE.vcf.gz && \\\n";
        $cmd.= "rm $mydir/*MuSE.txt $mydir/head.txt $mydir/txt.head.txt ";
        generateShell($shell,$cmd);
        $scripts{muse}{step2}{$sampid}=$shell;
        $snvvcf{$sampid}{muse}="$mydir/$sampid.MuSE.vcf.gz";
    }   
}

sub strelka{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/strelka";
        my $tbam = $hash{$sampid}{Tumor};
        my $nbam = $hash{$sampid}{Normal};
        my $shell="$outdir/$sampid/0.shell/a06.strelka/strelka.sh";
        my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        $cmd.= "$configs{strelkapl} --normal=$nbam --tumor=$tbam --ref=$configs{reference} --config=$configs{strelkapar} --output-dir=$mydir";
        generateShell($shell,$cmd);
        $scripts{strelka}{step1}{$sampid}=$shell;
        for my $mychr(@chrs){
            $shell = "$outdir/$sampid/0.shell/a06.strelka/runstrelka.$mychr.sh";
            # $cmd = "cd $mydir/$mychr\nmake";
            $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
            $cmd.= "if [ -d \"$mydir/$mychr\" ];then\n";
            $cmd.= "\tcd $mydir/$mychr && make\n";
            $cmd.= "fi";
            generateShell($shell,$cmd);
            $scripts{strelka}{step2}{$sampid}{$mychr}=$shell;
        }
        $shell="$outdir/$sampid/0.shell/a06.strelka/mergeStrelka.sh";
        $cmd = "less $mydir/chr1/results/passed.somatic.snvs.vcf | grep '#' > $mydir/head.snv.txt && \\\n";
        $cmd.= "less $mydir/chr1/results/passed.somatic.indels.vcf | grep '#' > $mydir/head.indel.txt && \\\n";
        $cmd.= "cat $mydir/*/results/passed.somatic.snvs.vcf | grep -v '#' | grep 'PASS' | cat $mydir/head.snv.txt - | gzip -f > $mydir/strelka.snv.vcf.gz && \\\n";
        $cmd.= "cat $mydir/*/results/passed.somatic.indels.vcf | grep -v '#' | grep 'PASS' | cat $mydir/head.indel.txt - | gzip -f > $mydir/$sampid.strelka.indel.vcf.gz && \\\n";
        $cmd.= "rm -r $mydir/chr* $mydir/head.snv.txt $mydir/head.indel.txt";
        generateShell($shell,$cmd);
        $scripts{strelka}{step3}{$sampid}=$shell;
        $snvvcf{$sampid}{strelka}="$mydir/strelka.snv.vcf.gz";
        $indelvcf{$sampid}{strelka}="$mydir/$sampid.strelka.indel.vcf.gz";
    }
}

sub mergeSnvInDel{
    my (%hash) = @_;
    for my $sampid(sort keys %snvvcf){
        my $mydir = "$outdir/$sampid/3.somatic/result";
        my $mutect = $snvvcf{$sampid}{mutect} if ($configs{mutect} =~ /true/i);
        my $muse = $snvvcf{$sampid}{muse} if ($configs{MuSE} =~ /true/i);
        my $mutect2 = $snvvcf{$sampid}{mutect2} if ($configs{mutect2} =~ /true/i);
        my $strelka = $snvvcf{$sampid}{strelka} if ($configs{strelka} =~ /true/i);
        my $strelka_indel = $indelvcf{$sampid}{strelka} if ($configs{strelka} =~ /true/i);
        my $strelka2 = $snvvcf{$sampid}{strelka2} if ($configs{strelka2} =~ /true/i);       
        my $strelka2_indel = $indelvcf{$sampid}{strelka2} if ($configs{strelka2} =~ /true/i);
        my $svaba = $indelvcf{$sampid}{svaba} if ($configs{svaba} =~ /true/i);
        for my $mychr(@chrs){
            system("mkdir -p $mydir/tmp") unless (-e "$mydir/tmp");
            my $shell = "$outdir/$sampid/0.shell/a10.mergeMut/$sampid.merge_$mychr.sh";
            my $cmd = "";
            $cmd.= "less $mutect | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.mutect.${mychr}.PASS.vcf && \\\n" if ($configs{mutect} =~ /true/i);
            $cmd.= "less $muse | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.muse.${mychr}.PASS.vcf && \\\n" if ($configs{MuSE} =~ /true/i);
            $cmd.= "less $mutect2 | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.mutect2.${mychr}.PASS.vcf && \\\n" if ($configs{mutect2} =~ /true/i);
            $cmd.= "less $strelka | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.strelka.${mychr}.PASS.vcf && \\\n" if ($configs{strelka} =~ /true/i);
            $cmd.= "less $strelka_indel | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.strelka_indel.${mychr}.PASS.vcf && \\\n" if ($configs{strelka} =~ /true/i);
            $cmd.= "less $strelka2 | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.strelka2.${mychr}.PASS.vcf && \\\n" if ($configs{strelka2} =~ /true/i);
            $cmd.= "less $strelka2_indel | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.strelka2_indel.${mychr}.PASS.vcf && \\\n" if ($configs{strelka2} =~ /true/i);
            $cmd.= "less $svaba | awk 'BEGIN{FS=\"\\t\"}(\$1~\"${mychr}_\" || \$1==\"$mychr\" || \$1~\"^#C\") && (\$7==\"PASS\" || \$7==\"FILTER\"){print \$0}' \> $mydir/tmp/${sampid}.svaba.${mychr}.PASS.vcf && \\\n" if ($configs{svaba} =~ /true/i);
            $cmd.= "$configs{mergeSnvpy}";
            $cmd.= " $mydir/tmp/${sampid}.mutect.${mychr}.PASS.vcf" if ($configs{mutect} =~ /true/i);
            $cmd.= " $mydir/tmp/${sampid}.mutect2.${mychr}.PASS.vcf" if ($configs{mutect2} =~ /true/i);
            $cmd.= " $mydir/tmp/${sampid}.strelka.${mychr}.PASS.vcf" if ($configs{strelka} =~ /true/i);
            $cmd.= " $mydir/tmp/${sampid}.strelka2.${mychr}.PASS.vcf" if ($configs{strelka2} =~ /true/i);
            $cmd.= " $mydir/tmp/${sampid}.muse.${mychr}.PASS.vcf" if ($configs{MuSE} =~ /true/i);
            $cmd.= " $mydir/tmp/${mychr}.merged_snv.unsorted.vcf $mydir/tmp/${mychr}.uniq_snv.unsorted.vcf && \\\n";
            $cmd.= "$configs{mergeIndelpy}";
            $cmd.= " $mydir/tmp/${sampid}.mutect2.${mychr}.PASS.vcf" if ($configs{mutect2} =~ /true/i);
            $cmd.= " $mydir/tmp/${sampid}.strelka_indel.${mychr}.PASS.vcf" if ($configs{strelka} =~ /true/i);
            $cmd.= " $mydir/tmp/${sampid}.strelka2_indel.${mychr}.PASS.vcf" if ($configs{strelka2} =~ /true/i);
            $cmd.= " $mydir/tmp/${sampid}.svaba.${mychr}.PASS.vcf" if ($configs{svaba} =~ /true/i);
            $cmd.= " $mydir/tmp/${mychr}.merged_indel.unsorted.vcf $mydir/tmp/${mychr}.uniq_indel.unsorted.vcf";
            generateShell($shell,$cmd);
            $scripts{mergeSnvInDel}{$sampid}{$mychr} = $shell;
        }
        my $shell = "$outdir/$sampid/0.shell/a10.mergeMut/$sampid.concat.sh"; 
        my $cmd = "less $mydir/tmp/chr*.merged_snv.unsorted.vcf > $mydir/merged_snv.unsorted.vcf && \\\n";
        $cmd.= "less $mydir/tmp/chr*.uniq_snv.unsorted.vcf > $mydir/uniq_snv.unsorted.vcf && \\\n";
        $cmd.= "less $mydir/tmp/chr*.merged_indel.unsorted.vcf > $mydir/merged_indel.unsorted.vcf && \\\n";
        $cmd.= "less $mydir/tmp/chr*.uniq_indel.unsorted.vcf > $mydir/uniq_indel.unsorted.vcf && \\\n";
        $cmd.= "$configs{vcfsort} $mydir/merged_snv.unsorted.vcf | gzip -c > $mydir/$sampid.RNA.somatic.merged_snv.vcf.gz && \\\n";
        $cmd.= "$configs{vcfsort} $mydir/uniq_snv.unsorted.vcf | gzip -c > $mydir/$sampid.somatic.RNA.uniq_snv.vcf.gz && \\\n";
        $cmd.= "$configs{vcfsort} $mydir/merged_indel.unsorted.vcf | gzip -c > $mydir/$sampid.RNA.somatic.merged_indel.vcf.gz && \\\n";
        $cmd.= "$configs{vcfsort} $mydir/uniq_indel.unsorted.vcf | gzip -c > $mydir/$sampid.somatic.RNA.uniq_indel.vcf.gz && \\\n";
        $cmd.= "rm -rf $mydir/tmp $mydir/*.unsorted.vcf";
        generateShell($shell,$cmd);
        $scripts{mergeSnvInDel}{$sampid}{concat} = $shell;

        $shell = "$outdir/$sampid/0.shell/a10.mergeMut/phase.sh";
        # $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        # $cmd.= "$configs{vep} --vcf -i $mydir/merged_snv.vcf -o $mydir/result.merged_snv.annot.input --offline --merged --cache --dir_cache $configs{vepdatabase} --assembly GRCh38 --dir_plugins $configs{vepplugins} --force_overwrite --use_given_ref $configs{reference} && \\\n";
        # $cmd.= "$configs{vep} --vcf -i $mydir/merged_indel.vcf -o $mydir/result.merged_indel.annot.input --offline --merged --cache --dir_cache $configs{vepdatabase} --assembly GRCh38 --dir_plugins $configs{vepplugins} --force_overwrite --use_given_ref $configs{reference}";
        $cmd = "perl $configs{ANNOVAR} $mydir/$sampid.RNA.somatic.merged_snv.vcf.gz $configs{ANNOVAR_refdir} -out $mydir/$sampid.RNA.somatic.merged_snv.annot $configs{annovar_par} && \\\n";
        $cmd.= "ls $mydir/$sampid.RNA.somatic.merged_snv.annot* | while read file; do $configs{bgzip} \$file; done && \\\n";
        $cmd.= "perl $configs{ANNOVAR} $mydir/$sampid.RNA.somatic.merged_indel.vcf.gz $configs{ANNOVAR_refdir} -out $mydir/$sampid.RNA.somatic.merged_indel.annot $configs{annovar_par} && \\\n";
        $cmd.= "ls $mydir/$sampid.RNA.somatic.merged_indel.annot* | while read file; do $configs{bgzip} \$file; done";
        generateShell($shell,$cmd);
        $scripts{phase}{$sampid} = $shell;
        $snvvcf{$sampid}{merge} = "$mydir/$sampid.RNA.somatic.merged_snv.vcf.gz";
        $indelvcf{$sampid}{merge} = "$mydir/$sampid.RNA.somatic.merged_indel.vcf.gz";
    }     
}

sub hlasomatic{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/hlasomatic";
        my $tbam = $hash{$sampid}{Tumor};
        my $nbam = $hash{$sampid}{Normal};
        my $shell1 = "$outdir/$sampid/0.shell/a11.hlasomatic/$sampid.normalHLA.sh";
        my $shell2 = "$outdir/$sampid/0.shell/a11.hlasomatic/$sampid.tumorHLA.sh";
        system("mkdir -p $mydir/Normal/Polysolver") unless (-e "$mydir/Normal/Polysolver");
        system("mkdir -p $mydir/Tumor/Polysolver") unless (-e "$mydir/Tumor/Polysolver");
        open OUT,"> $mydir/sort.bam.list";
        print OUT "SampleType\tPatientName\tNormal_bam\tTumor_bam\n";
        print OUT "SKCM\t$sampid\t$nbam\t$tbam\n";
        close OUT;
        my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        $cmd.= "source $configs{ploysolver_config} && \\\n";
        $cmd.= "$configs{polysolver} $nbam $configs{hlapar} $mydir/Normal/Polysolver && \\\n";
        $cmd.= "perl $configs{combinescript} $mydir/Normal/Polysolver $mydir/HLAdect.normal.csv && \\\n";
        $cmd.= "rm -rf $mydir/Normal/Polysolver/temp.*.bam";
        generateShell($shell1,$cmd);
        $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        $cmd.= "source $configs{ploysolver_config} && \\\n";
        $cmd.= "$configs{polysolver} $tbam $configs{hlapar} $mydir/Tumor/Polysolver && \\\n";
        $cmd.= "perl $configs{combinescript} $mydir/Tumor/Polysolver $mydir/HLAdect.tumor.csv && \\\n";
        $cmd.= "rm -rf $mydir/Tumor/Polysolver/temp.*.bam";
        generateShell($shell2,$cmd);
        $scripts{hlasomatic}{$sampid}{normal} = $shell1;
        $scripts{hlasomatic}{$sampid}{tumor} = $shell2;
        $hlatype{$sampid} = "$mydir/HLAdect.normal.csv";
    }

}

sub GATK_TumorOnly{
    my (%hash) = @_;
    for my $sampid(sort keys %hash){
        my $mydir = "$outdir/$sampid/3.somatic/tumor_only";
        system("mkdir -p $mydir/contamination/tmp") unless (-e "$mydir/contamination/tmp");
        system("mkdir -p $mydir/artifact/tmp") unless (-e "$mydir/artifact/tmp");
        system("mkdir -p $mydir/tmp") unless (-e "$mydir/tmp");
        my $tbam = $hash{$sampid}{Tumor};
        ## Contamination Estimate
        my $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.contamination.sh";
        my $cmd = "#$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\" GetPileupSummaries -I $tbam -V $configs{GetPileupSummaries} -L $configs{reference_bed} -O $mydir/contamination/$sampid\_Tumor.pileups.table && \\\n";
        $cmd.= "#$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\" CalculateContamination -I $mydir/contamination/$sampid\_Tumor.pileups.table -O $mydir/contamination/$sampid.contamination.table";
        generateShell($shell,$cmd);
        $scripts{tumor_only}{$sampid}{contamination} = $shell;
        ## Artifact calculation
        # $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.artifact.sh";
        # $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        # $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/artifact/tmp\" CollectSequencingArtifactMetrics -R $configs{reference} -I $tbam --FILE_EXTENSION \".txt\" -O $mydir/artifact/$sampid.artifact";
        # generateShell($shell,$cmd);
        # $scripts{tumor_only}{$sampid}{artifact} = $shell;
        ## Mutect2
        my @vcfs;
        for my $mychr(@chrs){
            my $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.call.$mychr.sh";
            my $thistmp = "$mydir/$mychr.tmp";
            system("mkdir -p $thistmp") unless (-e "$thistmp");
            $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
            $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$thistmp\" Mutect2 -R $configs{reference} -I $tbam -tumor $sampid\_Tumor";
            $cmd.= " --germline-resource $configs{gnomAD}" if($configs{gnomAD});
            $cmd.= " --panel-of-normals $configs{panel_of_normal}" if($configs{panel_of_normal});
            $cmd.= " -L $configs{reference_block_bed_path}/$mychr.bed -O $mydir/$sampid.mutect2.$mychr.vcf.gz && \\\n";
            $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterMutectCalls -V $mydir/$sampid.mutect2.$mychr.vcf.gz -O $mydir/$sampid.mutect2.$mychr.filter1.vcf.gz -R $configs{reference} && \\\n";
            $cmd.= "rm $mydir/$sampid.mutect2.$mychr.vcf.gz*";
            generateShell($shell,$cmd);
            $scripts{tumor_only}{step1}{$sampid}{$mychr} = $shell;
            push @vcfs, "$mydir/$sampid.mutect2.$mychr.filter1.vcf.gz ";
        }
        $shell = "$outdir/$sampid/0.shell/a05.mutect2/mergevcf.mutect2.sh";
        $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
        $cmd.= "$configs{vcfconcat} @vcfs | $configs{vcfsort} -t $mydir > $mydir/$sampid.mutect2.filter1.vcf && \\\n";
        $cmd.= "$configs{bgzip} $mydir/$sampid.mutect2.filter1.vcf && \\\n";
        $cmd.= "$configs{tabix} $mydir/$sampid.mutect2.filter1.vcf.gz && \\\n";
        $cmd.= "rm @vcfs";
        generateShell($shell,$cmd);
        $scripts{tumor_only}{step2}{$sampid} = $shell;
        ## Filter Mutect Calls
        $shell = "$outdir/$sampid/0.shell/a05.mutect2/mutect2.bqsr.filter.sh";
        $cmd = "#$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterMutectCalls -V $mydir/$sampid.mutect2.vcf.gz -contamination-table $mydir/contamination/$sampid.contamination.table -O $mydir/$sampid.mutect2.filter1.vcf.gz -R $configs{reference} && \\\n";
        $cmd.= "gzip -dc $mydir/$sampid.mutect2.filter1.vcf.gz | awk 'BEGIN{FS=\"\\t\"}(\$1~\"^#\" || \$7==\"PASS\"){print \$0}' | gzip -c \> $mydir/$sampid.tumor_only.somatic.vcf.gz && \\\n";
        $cmd.= "perl $configs{ANNOVAR} $mydir/$sampid.tumor_only.somatic.vcf.gz $configs{ANNOVAR_refdir} -out $mydir/$sampid.tumor_only.annot $configs{annovar_par} && \\\n";
        $cmd.= "ls $mydir/$sampid.tumor_only.annot* | while read file; do $configs{bgzip} \$file; done && \\\n";
        $cmd.= "rm -rf $mydir/chr*.tmp $mydir/$sampid.mutect2.chr* ";
        generateShell($shell,$cmd);
        $scripts{tumor_only}{$sampid}{filter} = $shell;
        $snvvcf{$sampid}{merge} = "$mydir/$sampid.tumor_only.somatic.vcf.gz";
        $indelvcf{$sampid}{merge} = "$mydir/$sampid.tumor_only.somatic.vcf.gz";
    }
}

sub kraken2{
    my (%hash) = @_;
    for my $sampid(keys %hash){
        for my $dtype(keys %{$hash{$sampid}}){
            my $shell = "$outdir/$sampid/0.shell/a12.kraken2/a12.kraken2_$dtype.sh";
            my $fastq1 = $hash{$sampid}{$dtype}{fq1};
            my $fastq2 = $hash{$sampid}{$dtype}{fq2};
            my $cmd = "mkdir -p $outdir/$sampid/4.microbiome/$dtype && \\\n";
            $cmd.= "$configs{kraken2} --db $configs{Kraken2DB} --threads $configs{Kraken2_threads} $configs{Kraken2par}";
            $cmd.= " --paired --classified-out $outdir/$sampid/4.microbiome/$dtype/cseq#.fq";
            $cmd.= " --output $outdir/$sampid/4.microbiome/$dtype/${sampid}_${dtype}_kraken2.out";
            $cmd.= " --report $outdir/$sampid/4.microbiome/$dtype/${sampid}_${dtype}_kraken2.report $fastq1 $fastq2";
            generateShell($shell, $cmd);
            $scripts{kraken2}{$sampid}{$dtype} = $shell;
        }
    }
}

sub edgeList{
    my ($day) = @_;
    my $mydir = "$outdir/shell_run";
    open EDGE, ">$outdir/shell_run/edge.$configs{projectname}.$day.list" or die $!;
    %fqinputs = %bqsrbam if ($configs{alignment} =~ /false/i);
    for my $sampid(sort keys %fqinputs){
        for my $dtype(sort keys %{$fqinputs{$sampid}}){
            # clean fastq
            if ($configs{fastqclean} =~ /true/i){        
                for my $lane(sort keys %{$scripts{clean}{$sampid}{$dtype}}){
                    for my $num(sort keys %{$scripts{clean}{$sampid}{$dtype}{$lane}}){
                        if($configs{split_num} > 1 && $splitfqs{$sampid}{$dtype}{$lane}{fqstat_exist} =~ /true/i){
                            ## split fq
                            my $t = $num % $configs{split_num};
                            my $zui = ($num - $t) / $configs{split_num}; 
                            print EDGE "$scripts{split}{$sampid}{$dtype}{$lane}{$zui.1}:1G:1CPU $scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU\n";
                            print EDGE "$scripts{split}{$sampid}{$dtype}{$lane}{$zui.2}:1G:1CPU $scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU\n";   
                        }
                        if($configs{fastqcorrect} =~ /true/i){      
                            for my $i(sort keys %{$scripts{correct}{$sampid}{$dtype}}){
                                print EDGE "$scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU $scripts{correct}{$sampid}{$dtype}{$i}:20G:4CPU\n";
                            }
                        }else{
                            print EDGE "$scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU $scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q}\n" if($configs{alignment} =~ /true/i || $configs{addreadgroup} =~ /true/i);
                        }
                    }
               }
            }
            # alignment
            if($configs{fastqcorrect} =~ /true/i){
                for my $i(sort keys %{$scripts{correct}{$sampid}{$dtype}}){
                    print EDGE "$scripts{correct}{$sampid}{$dtype}{$i}:20G:4CPU $scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q}\n" if($configs{alignment} =~ /true/i || $configs{addreadgroup} =~ /true/i);
                }
            }
            if($configs{bampostprocess} =~ /true/i){
                for my $mychr(@chrs){
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{bampostprocess}{$sampid}{$dtype}{$mychr}:10G:2CPU\n";
                    print EDGE "$scripts{bampostprocess}{$sampid}{$dtype}{$mychr}:10G:2CPU $scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU\n";
                    print EDGE "$scripts{bampostprocess}{$sampid}{$dtype}{$mychr}:10G:2CPU $scripts{GatherBQSR}{$sampid}{$dtype}:5G:1CPU\n";
                }
            }
            # BAM QC
            if($configs{alignment} =~ /true/i || $configs{addreadgroup} =~ /true/i){
                if($configs{bampostprocess} =~ /false/i){
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{coverage}{$sampid}{$dtype}:10G:1CPU\n" if($configs{bamqc} =~ /true/i);
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{BAMqc}{$sampid}:8G:1CPU\n" if($configs{bamqc} =~ /true/i && $configs{qcvcf});
                }else{
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{coverage}{$sampid}{$dtype}:10G:1CPU\n" if($configs{bamqc} =~ /true/i);
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{BAMqc}{$sampid}:8G:1CPU\n" if($configs{bamqc} =~ /true/i && $configs{qcvcf});
                }
            }
            for my $mychr(@chrs){
                ## Mutect
                if($configs{mutect} =~ /true/i){
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{mutect}{step1}{$sampid}{$mychr}:10G:1CPU\n" if($configs{bampostprocess} =~ /false/i);
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{mutect}{step1}{$sampid}{$mychr}:10G:1CPU\n" if($configs{bampostprocess} =~ /true/i);
                    print EDGE "$scripts{mutect}{step1}{$sampid}{$mychr}:10G:1CPU $scripts{mutect}{step2}{$sampid}:1G:1CPU\n" if($dtype =~ /Normal/i);
                }
                # MuSe
                if($configs{MuSE} =~ /true/i){
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{MuSE}{step1}{$sampid}{$mychr}:8G:1CPU\n" if($configs{bampostprocess} =~ /false/i);
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{MuSE}{step1}{$sampid}{$mychr}:8G:1CPU\n" if($configs{bampostprocess} =~ /true/i);
                    print EDGE "$scripts{MuSE}{step1}{$sampid}{$mychr}:8G:1CPU $scripts{muse}{step2}{$sampid}:20G:1CPU\n" if ($dtype =~ /Normal/i);
                }
                # SvaBA
                if($configs{svaba} =~ /true/i){
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{svaba}{$sampid}{$mychr}:10G:1CPU\n" if($configs{bampostprocess} =~ /false/i);
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{svaba}{$sampid}{$mychr}:10G:1CPU\n" if($configs{bampostprocess} =~ /true/i);
                    print EDGE "$scripts{svaba}{$sampid}{$mychr}:10G:1CPU $scripts{svaba}{$sampid}{merge}:1G:1CPU\n" if($dtype =~ /Normal/i);
                }
                # Strelka
                if($configs{strelka} =~ /true/i && $dtype =~ /Normal/i){
                    print EDGE "$scripts{strelka}{step1}{$sampid}:8G:1CPU $scripts{strelka}{step2}{$sampid}{$mychr}:8G:1CPU\n";
                    print EDGE "$scripts{strelka}{step2}{$sampid}{$mychr}:8G:1CPU $scripts{strelka}{step3}{$sampid}:1G:1CPU\n";
                }
                # Mutect2
                if($configs{mutect2} =~ /true/i){
                    # print EDGE "$scripts{mutect2}{$sampid}{PoN}:15G:2CPU $scripts{mutect2}{step1}{$sampid}{$mychr}:20G:2CPU\n" unless($configs{panel_of_normal});
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{mutect2}{step1}{$sampid}{$mychr}:20G:2CPU\n" if((($configs{alignment} =~ /true/i || $configs{addreadgroup} =~ /true/i) && $configs{bampostprocess} =~ /false/i) && $configs{panel_of_normal});
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{mutect2}{step1}{$sampid}{$mychr}:20G:2CPU\n"  if($configs{bampostprocess} =~ /true/i);
                    print EDGE "$scripts{mutect2}{step1}{$sampid}{$mychr}:20G:2CPU $scripts{mutect2}{step2}{$sampid}:5G:1CPU\n" if($dtype =~ /Normal/);
                }
                # tumor only mode
                if($configs{tumor_only} =~ /true/i && $dtype =~ /Tumor/){
                    # print EDGE "$scripts{tumor_only}{$sampid}{PoN}:15G:2CPU $scripts{tumor_only}{step1}{$sampid}{$mychr}:8G:1CPU\n" unless($configs{panel_of_normal});
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{tumor_only}{step1}{$sampid}{$mychr}:8G:1CPU\n" if($configs{bampostprocess} =~ /false/i);
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{tumor_only}{step1}{$sampid}{$mychr}:8G:1CPU\n" if($configs{bampostprocess} =~ /true/i);
                    print EDGE "$scripts{tumor_only}{step1}{$sampid}{$mychr}:8G:1CPU $scripts{tumor_only}{step2}{$sampid}:5G:1CPU\n";
                }
            }

            ## strelka (BQSR bam)
            if($configs{strelka} =~ /true/i){
                print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{strelka}{step1}{$sampid}:8G:1CPU\n" if($configs{bampostprocess} =~ /false/i);
                print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{strelka}{step1}{$sampid}:8G:1CPU\n" if($configs{bampostprocess} =~ /true/i);
            }
            ## Strelka2 (BQSR bam)
            if($configs{strelka2} =~ /true/i){
                print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{strelka2}{$sampid}:6G:4CPU\n" if($configs{bampostprocess} =~ /false/i);
                print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{strelka2}{$sampid}:6G:4CPU\n" if($configs{bampostprocess} =~ /true/i);
            }
            ## Mutect2 (BQSR bam)
            if($configs{mutect2} =~ /true/i){
                # if($configs{bampostprocess} =~ /false/i){
                    print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{mutect2}{$sampid}{contamination}:1G:1CPU\n" if($configs{bampostprocess} =~ /true/i);
                    print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{mutect2}{$sampid}{contamination}:1G:1CPU\n" if($configs{bampostprocess} =~ /false/i);
                    # print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{mutect2}{$sampid}{artifact}:4G:1CPU\n" if($dtype =~ /Tumor/i);
                    # print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{mutect2}{$sampid}{PoN}:15G:2CPU\n" unless($configs{panel_of_normal});
                # }
                if($dtype =~ /Tumor/i){
                    print EDGE "$scripts{mutect2}{$sampid}{contamination}:1G:1CPU $scripts{mutect2}{$sampid}{filter}:1G:1CPU\n";
                    print EDGE "$scripts{mutect2}{step2}{$sampid}:5G:1CPU $scripts{mutect2}{$sampid}{filter}:1G:1CPU\n";
                    # print EDGE "$scripts{mutect2}{$sampid}{artifact}:4G:1CPU $scripts{mutect2}{$sampid}{filter}:1G:1CPU\n";
                }
            }
            ## tumor only (BQSR bam)
            if($configs{tumor_only} =~ /true/i){
                print EDGE "$scripts{BamMerge}{$sampid}{$dtype}:10G:2CPU $scripts{tumor_only}{$sampid}{contamination}:1G:1CPU\n" if($configs{bampostprocess} =~ /true/i);
                print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{tumor_only}{$sampid}{contamination}:1G:1CPU\n" if($configs{bampostprocess} =~ /false/i);
                if($dtype =~ /Tumor/i){
                    print EDGE "$scripts{tumor_only}{$sampid}{contamination}:1G:1CPU $scripts{tumor_only}{$sampid}{filter}:35G:1CPU\n";
                    print EDGE "$scripts{tumor_only}{step2}{$sampid}:5G:1CPU $scripts{tumor_only}{$sampid}{filter}:35G:1CPU\n";
                    # print EDGE "$scripts{tumor_only}{$sampid}{artifact}:4G:1CPU $scripts{tumor_only}{$sampid}{filter}:35G:1CPU\n";
                }
            }
            if ($configs{somatic} =~ /true/i && $dtype =~ /Normal/i){
                for my $mychr(@chrs){
                    print EDGE "$scripts{mutect}{step2}{$sampid}:1G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{mutect} =~ /true/i);
                    print EDGE "$scripts{mutect2}{$sampid}{filter}:1G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n";
                    print EDGE "$scripts{strelka}{step3}{$sampid}:6G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{strelka} =~ /true/i);
                    print EDGE "$scripts{muse}{step2}{$sampid}:20G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{MuSE} =~ /true/i);
                    print EDGE "$scripts{svaba}{$sampid}{merge}:1G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{svaba} =~/true/i);
                    print EDGE "$scripts{strelka2}{$sampid}:6G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{strelka2}=~/true/i);
                    print EDGE "$scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU $scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU\n";
                }
                print EDGE "$scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU $scripts{phase}{$sampid}:2G:2CPU\n";
            }
            ##hla detecion
            unless(-e $configs{hlalist}){
                print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{hlasomatic}{$sampid}{normal}:4G:1CPU\n" if($configs{hlasomatic} =~ /true/i);
                print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{hlasomatic}{$sampid}{tumor}:4G:1CPU\n" if($configs{hlasomatic} =~ /true/i);
            }
            ## microbiome detection
            print EDGE "$scripts{alignment}{$sampid}{$dtype}:$configs{alignment_q} $scripts{kraken2}{$sampid}{$dtype}:$configs{kraken2_q}\n" if($configs{kraken} =~ /true/i && $configs{alignment} =~ /true/i);
        }
    }
    close EDGE;
}
