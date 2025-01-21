#!/usr/bin/perl -w
use strict;
use POSIX;
use File::Basename;
use Data::Dumper;

@ARGV == 3 or die "Usage: perl $0 <Fastq.list> <Configs> <Outdir>\n";
my ($fqlist,$config,$outdir) = @ARGV;
system("mkdir -p $outdir") unless (-e "$outdir");
#######################################################################
my (%configs,%scripts);
my (%fqinputs,%cleanfqs,%splitfqs,%cleanfastq,%phred,%hlatype,%unmapfqs);
my (%bam,%realnbam,%realnbam2,%bqsrbam,%rnabam);
my (%snvvcf,%indelvcf,%fusionresult,%exprresult);

my @MegaBOLThead = qw{SampleName Read1 Read2};
my @chrs;
########################### Main Parameters ###########################

&readConfigs($config);
if($configs{tumor_only} =~ /true/i){
	$configs{svaba} = "false";
	$configs{MuSE} = "false";
	$configs{strelka} = "false";
	$configs{strelka2} = "false";
	$configs{mutect} = "false";
	$configs{mutect2} = "false";
	$configs{varscan2} = "false";
	$configs{FACTERA} = "false";
	$configs{manta} = "false";
	$configs{gatkcnv} = "false";
	$configs{FACETS} = "false";
	$configs{msi} = "false";
}

open LIST,"<$configs{reference_block}" or die "$!";
while(my $line = <LIST>){
	chomp $line;
	push @chrs, $line;
}
close LIST;
# my @chrs = (1..22,"X","Y");
&readFqlist($fqlist);
&makeDir(\%fqinputs) if($configs{alignment} =~ /true/i);
&makeDir(\%bqsrbam) if($configs{alignment} =~ /false/i);
&runsplit(\%fqinputs) if ($configs{alignment} =~ /true/i and $configs{split_num} > 1);

if ($configs{fastqclean} =~ /true/i){
	if($configs{split_num} > 1){&FastqClean(\%splitfqs);}else{&FastqClean(\%fqinputs);}
}else{
	&FastqLink(\%fqinputs);
}

##allign BQSR HaplotypeCaller
&MegaBOLT(%cleanfqs) if($configs{alignment} =~ /true/i);

##QC
&megaboltqc(%bqsrbam) if($configs{megaboltqc} =~ /true/i);

# metaphlan
&metaphlan(%unmapfqs) if($configs{metaphlan} =~ /true/i);
# Kraken2
&kraken2(%unmapfqs) if($configs{kraken} =~ /true/i);

##call SV
&Manta(%bqsrbam) if($configs{manta} =~ /true/i);
&svaba(%bqsrbam) if($configs{svaba} =~ /true/i);

## somatic mutation calling
&mutectSNV(%bqsrbam) if ($configs{mutect} =~ /true/i);
&mutect2(%bqsrbam) if ($configs{mutect2} =~ /true/i);
&strelka(%bqsrbam) if ($configs{strelka} =~ /true/i);
&strelka2(%bqsrbam) if ($configs{strelka2} =~ /true/i);
&varscan2(%bqsrbam) if ($configs{varscan2} =~ /true/i);
&MuSE(%bqsrbam) if ($configs{MuSE} =~ /true/i);
&mergeSnvInDel() if ($configs{mutect} =~ /true/i and $configs{mutect2} =~ /true/i);

## call SCNV
&FACETS(%bqsrbam) if ($configs{FACETS} =~ /true/i);
&gatkcnv1(%bqsrbam) if($configs{gatkcnv} =~ /true/i);

## gene fusion calling
&integrate(\%rnabam,\%bqsrbam) if ($configs{integrate} =~ /true/i);
&starfusion() if ($configs{starfusion} =~ /true/i);

## expression calculation
&rsemstar(\%bqsrbam) if ($configs{rsemstar} =~ /true/i);

## MSI
&msi(%bqsrbam) if ($configs{msi} =~ /true/i);

#HLA type
#&HLAtyping()
if(-e $configs{hlalist}){
	&readHLAlist();
}
else{
	&hlasomatic(%bqsrbam) if($configs{hlasomatic} =~ /true/i);
	&HLAminer(%cleanfqs) if($configs{hlaminer} =~ /true/i);
}

## medicine annotation
&targetDrug() if ($configs{targetDrug} =~ /true/i);
&resistantDrug() if ($configs{resistant} =~ /true/i);
&chemicalDrug(%bqsrbam) if ($configs{chemicalDrug} =~ /true/i && $configs{alignment} =~ /false/i); # if input is bam not fastq ###

## immupathway
&immupathway() if ($configs{immupathway} =~ /true/i);

# Generate day info
chomp(my $day = `date +%Y%m%d`);
## run scripts
&edgeList($day);
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

sub readFqlist{
	my ($inputs) = @_;
	open INPUT, "<$inputs" or die "$!";
	my $num = 0;
	while(my $line = <INPUT>){
		chomp $line;
		my ($sampid,$dtype,$phred,$lane,$reads) = split /\s+/,$line;
		$dtype = "Tumor" if($dtype =~ /Tumor/i);
		$dtype = "Normal" if($dtype =~ /Normal/i);
		if ($configs{alignment} =~ /false/i){
			$bqsrbam{$sampid}{$dtype} = $reads; ## bam file input
		}else{
			my ($fq1,$fq2) = split /\,/, $reads;
			push @{$fqinputs{$sampid}{$dtype}{$lane}{fq1}}, $fq1;
			push @{$fqinputs{$sampid}{$dtype}{$lane}{fq2}}, $fq2;
		}
		$phred{$sampid}{$dtype} = $phred;
	}
	close INPUT;
}

sub makeDir{
	my ($hash) = @_;
	my %hash = %{$hash};
	system("mkdir -p $outdir/shell_run") unless (-e "$outdir/shell_run");
	for my $sampid(keys %hash){
		system("mkdir -p $outdir/$sampid/0.shell") unless (-e "$outdir/$sampid/0.shell");
		system("mkdir -p $outdir/$sampid/1.clean") unless (-e "$outdir/$sampid/1.clean");
		system("mkdir -p $outdir/$sampid/2.megaBOLT") unless (-e "$outdir/$sampid/2.megaBOLT");
		system("mkdir -p $outdir/$sampid/3.SV") unless (-e "$outdir/$sampid/3.SV");
		system("mkdir -p $outdir/$sampid/4.SCNV") unless (-e "$outdir/$sampid/4.SCNV");
		system("mkdir -p $outdir/$sampid/5.somatic") unless (-e "$outdir/$sampid/5.somatic");
		system("mkdir -p $outdir/$sampid/6.fusion") unless (-e "$outdir/$sampid/6.fusion");
		system("mkdir -p $outdir/$sampid/7.expression") unless (-e "$outdir/$sampid/7.expression");
		system("mkdir -p $outdir/$sampid/8.neoantigen") unless (-e "$outdir/$sampid/8.neoantigen");
		system("mkdir -p $outdir/$sampid/9.medicine") unless (-e "$outdir/$sampid/9.medicine");
		system("mkdir -p $outdir/$sampid/10.immupathway") unless (-e "$outdir/$sampid/10.immupathway");
		system("mkdir -p $outdir/$sampid/11.microbiome") unless (-e "mkdir -p $outdir/$sampid/11.microbiome");
		## shell
		system("mkdir -p $outdir/$sampid/0.shell/a01.clean") unless (-e "$outdir/$sampid/0.shell/a01.clean");
		system("mkdir -p $outdir/$sampid/0.shell/a02.megaBOLT") unless (-e "$outdir/$sampid/0.shell/a02.megaBOLT");
		system("mkdir -p $outdir/$sampid/0.shell/a05.mutect") unless (-e "$outdir/$sampid/0.shell/a05.mutect");
		system("mkdir -p $outdir/$sampid/0.shell/a06.mutect2") unless (-e "$outdir/$sampid/0.shell/a06.mutect2");
		system("mkdir -p $outdir/$sampid/0.shell/a07.strelka") unless (-e "$outdir/$sampid/0.shell/a07.strelka");
		system("mkdir -p $outdir/$sampid/0.shell/a08.strelka2") unless (-e "$outdir/$sampid/0.shell/a08.strelka2");
		system("mkdir -p $outdir/$sampid/0.shell/a09.FACETS") unless (-e "$outdir/$sampid/0.shell/a09.FACETS");
		system("mkdir -p $outdir/$sampid/0.shell/a10.integrate") unless (-e "$outdir/$sampid/0.shell/a10.integrate");
		system("mkdir -p $outdir/$sampid/0.shell/a11.starfusion") unless (-e "$outdir/$sampid/0.shell/a11.starfusion");
		system("mkdir -p $outdir/$sampid/0.shell/a12.svaba") unless (-e "$outdir/$sampid/0.shell/a12.svaba");
		system("mkdir -p $outdir/$sampid/0.shell/a13.rsem") unless (-e "$outdir/$sampid/0.shell/a13.rsem");
		system("mkdir -p $outdir/$sampid/0.shell/a14.mergeMut") unless (-e "$outdir/$sampid/0.shell/a14.mergeMut");
		system("mkdir -p $outdir/$sampid/0.shell/a15.neoantigen") unless (-e "$outdir/$sampid/0.shell/a15.neoantigen");
		system("mkdir -p $outdir/$sampid/0.shell/a16.medicine") unless (-e "$outdir/$sampid/0.shell/a16.medicine");
		system("mkdir -p $outdir/$sampid/0.shell/a17.immupathway") unless (-e "$outdir/$sampid/0.shell/a17.immupathway");
		system("mkdir -p $outdir/$sampid/0.shell/a18.hlasomatic") unless (-e "$outdir/$sampid/0.shell/a18.hlasomatic");
		system("mkdir -p $outdir/$sampid/0.shell/a19.gatkcnv") unless (-e "$outdir/$sampid/0.shell/a19.gatkcnv");
		system("mkdir -p $outdir/$sampid/0.shell/a20.Manta") unless (-e "$outdir/$sampid/0.shell/a20.Manta");
		system("mkdir -p $outdir/$sampid/0.shell/a21.msi") unless (-e "$outdir/$sampid/0.shell/a21.msi");
		system("mkdir -p $outdir/$sampid/0.shell/a22.varscan2") unless (-e "$outdir/$sampid/0.shell/a22.varscan2");
		system("mkdir -p $outdir/$sampid/0.shell/a23.MuSE") unless (-e "$outdir/$sampid/0.shell/a23.MuSE");
		system("mkdir -p $outdir/$sampid/0.shell/a24.metaphlan") unless (-e "$outdir/$sampid/0.shell/a24.metaphlan");
		system("mkdir -p $outdir/$sampid/0.shell/a25.kraken2")unless(-e "$outdir/$sampid/0.shell/a25.kraken2");
		## somatic mutation
		system("mkdir -p $outdir/$sampid/5.somatic/mutect") unless (-e "$outdir/$sampid/5.somatic/mutect");
		system("mkdir -p $outdir/$sampid/5.somatic/mutect2") unless (-e "$outdir/$sampid/5.somatic/mutect2");
		system("mkdir -p $outdir/$sampid/5.somatic/strelka") unless (-e "$outdir/$sampid/5.somatic/strelka");
		system("mkdir -p $outdir/$sampid/5.somatic/strelka2") unless (-e "$outdir/$sampid/5.somatic/strelka2");
		system("mkdir -p $outdir/$sampid/5.somatic/svaba") unless (-e "$outdir/$sampid/5.somatic/svaba");
		system("mkdir -p $outdir/$sampid/5.somatic/result") unless (-e "$outdir/$sampid/5.somatic/result");
		system("mkdir -p $outdir/$sampid/5.somatic/hlasomatic") unless (-e "$outdir/$sampid/5.somatic/hlasomatic");
		system("mkdir -p $outdir/$sampid/5.somatic/msi") unless (-e "$outdir/$sampid/5.somatic/msi");
		system("mkdir -p $outdir/$sampid/5.somatic/varscan2") unless (-e "$outdir/$sampid/5.somatic/varscan2");
		system("mkdir -p $outdir/$sampid/5.somatic/MuSE") unless (-e "$outdir/$sampid/5.somatic/MuSE");
		## gene fusion
		system("mkdir -p $outdir/$sampid/6.fusion/integrate") unless (-e "$outdir/$sampid/6.fusion/integrate");
		system("mkdir -p $outdir/$sampid/6.fusion/starfusion") unless (-e "$outdir/$sampid/6.fusion/starfusion");
		system("mkdir -p $outdir/$sampid/6.fusion/edicofusion") unless (-e "$outdir/$sampid/6.fusion/edicofusion");
		## gene expression
		system("mkdir -p $outdir/$sampid/7.expression/RSEMstar") unless (-e "$outdir/$sampid/7.expression/RSEMstar");
		##cnv
		system("mkdir -p $outdir/$sampid/4.SCNV/FACETS") unless (-e "$outdir/$sampid/4.SCNV/FACETS");
		system("mkdir -p $outdir/$sampid/4.SCNV/gatkcnv") unless (-e "$outdir/$sampid/4.SCNV/gatkcnv");
		## medicine
		system("mkdir -p $outdir/$sampid/9.medicine/chemicalDrug") unless (-e "$outdir/$sampid/9.medicine/chemicalDrug");
		system("mkdir -p $outdir/$sampid/9.medicine/targetDrug") unless (-e "$outdir/$sampid/9.medicine/targetDrug");
		system("mkdir -p $outdir/$sampid/9.medicine/resistantDrug") unless (-e "$outdir/$sampid/9.medicine/resistantDrug");
		##qc
		system("mkdir -p $outdir/$sampid/2.megaBOLT/QC") unless (-e "$outdir/$sampid/2.megaBOLT/QC");
		
		
		if($configs{alignment} =~ /true/i){
			for my $dtype(keys %{$hash{$sampid}}){
				for my $lane(keys %{$hash{$sampid}{$dtype}}){
					system("mkdir -p $outdir/$sampid/1.clean/$dtype/$lane") unless (-e "$outdir/$sampid/1.clean/$dtype/$lane");
				}
			}
		}
		## Target region
		# `ln -s $configs{USER_TR} $outdir/$sampid/TR` unless (-e "$outdir/$sampid/TR");
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
	print OUT "$content && \\\n" if($content);
	print OUT "echo ==========end at : `date` ========== && \\\n";
	print OUT "echo $finish_string 1>&2 && \\\n";
	print OUT "echo $finish_string > $output_shell.sign\n";
	print OUT "qstat -j $shell_name > $output_shell.log\n";
	close OUT;
}

sub runsplit{
	my ($hash) = @_;
	my %hash = %{$hash};
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
	my ($hash) = @_;
	my %hash = %{$hash};
	for my $sampid(sort keys %hash){
		for my $dtype(sort keys %{$hash{$sampid}}){
			if ($dtype =~ /Tumor/i or $dtype =~ /Normal/i){
				$cleanfqs{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$dtype.csv";
				open ELIST, ">$outdir/$sampid/2.megaBOLT/$dtype.csv" or die "$!";
			}
			for my $lane(sort keys %{$hash{$sampid}{$dtype}}){
				my @fq1 = @{$hash{$sampid}{$dtype}{$lane}{fq1}};
				my @fq2 = @{$hash{$sampid}{$dtype}{$lane}{fq2}};
				for my $i(0..$#fq1){
					my ($out_fq1, $out_fq2);
					$out_fq1 = "$outdir/$sampid/1.clean/$dtype/$lane/fq$i/clean_1.$i.fq.gz";
					$out_fq2 = "$outdir/$sampid/1.clean/$dtype/$lane/fq$i/clean_2.$i.fq.gz";
					system("mkdir -p $outdir/$sampid/1.clean/$dtype/$lane/fq$i") unless (-e "$outdir/$sampid/1.clean/$dtype/$lane/fq$i");
					my $shell = "$outdir/$sampid/0.shell/a01.clean/a01.cleanfq_$sampid\_$dtype\_$lane\_$i.sh";
					my $cmd;
					if($configs{soapnukeclean} =~ /false/){
						$cmd = "$configs{fastp} -i $fq1[$i] -o $out_fq1";
						$cmd.= " -I $fq2[$i] -O $out_fq2" if($fq2[$i]);
						$cmd.= " -j $outdir/$sampid/1.clean/$dtype/$lane/fq$i/$sampid.$dtype.$lane.fastp.json";
						$cmd.= " -h $outdir/$sampid/1.clean/$dtype/$lane/fq$i/$sampid.$dtype.$lane.fastp.html";
						$cmd.= " $configs{fastpparameter}";
					}else{
						$cmd = "$configs{soapnuketool} filter -1 $fq1[$i] -C clean_1.$i.fq.gz";
						$cmd.= " -2 $fq2[$i] -D clean_2.$i.fq.gz" if($fq2[$i]);
						$cmd.= " $configs{soapnukeparameter} -o $outdir/$sampid/1.clean/$dtype/$lane/fq$i";
					}
					if($configs{split_num} > 1 && $hash{$sampid}{$dtype}{$lane}{fqstat_exist} =~ /true/i){
						$cmd.= "&& \\\nrm -rf $fq1[$i]";
						$cmd.= " && \\\nrm -rf $fq2[$i]" if($fq2[$i]);

					}
					generateShell($shell,$cmd);
					$scripts{clean}{$sampid}{$dtype}{$lane}{$i}=$shell;
					push @{$cleanfastq{$sampid}{$dtype}{fq1}}, $out_fq1;
					push @{$cleanfastq{$sampid}{$dtype}{fq2}}, $out_fq2;
				}
			}
			my $mega_fq1 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq1}});
			my $mega_fq2 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq2}});
			print ELIST join("\t","$sampid\-$dtype",$mega_fq1);
			print ELIST "\t$mega_fq2" unless($mega_fq2 =~ /,{2,}/);
			print ELIST "\n";
			close ELIST;
		}
	}
}

sub FastqLink{
	my ($hash) = @_;
	my %hash = %{$hash};
	for my $sampid(sort keys %hash){
		for my $dtype(sort keys %{$hash{$sampid}}){
			if ($dtype =~ /Tumor/i or $dtype =~ /Normal/i){
				$cleanfqs{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$dtype.csv";
				open ELIST, ">$outdir/$sampid/2.megaBOLT/$dtype.csv" or die "$!";
			}
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
			my $mega_fq1 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq1}});
			my $mega_fq2 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq2}});
			print ELIST join("\t","$sampid\-$dtype",$mega_fq1);
			print ELIST "\t$mega_fq2" unless($mega_fq2 =~ /,{2,}/);
			print ELIST "\n";
			close ELIST;
		}
	}
}

sub MegaBOLT{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			my $Suffix = "mm2";
			$Suffix = "bwa" if ($configs{BWA} =~ /true/i);
			$bqsrbam{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.bam";
			$unmapfqs{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$sampid-$dtype/unmap.fq.gz";
			my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/a02.megaBOLT_$dtype.sh";
			# my $cmd = "MegaBOLT --type basic --runtype WGS --outputdir $outdir/$sampid/2.megaBOLT/ --bwa 1";
			my $cmd = "MegaBOLT --type basic --runtype WGS --outputdir $outdir/$sampid/2.megaBOLT/";
			$cmd.= " --bwa 1" if ($configs{BWA} =~ /true/i); 
			# $cmd.= " --list $outdir/$sampid/2.megaBOLT/Tumor/Tumor.csv --list2 $outdir/$sampid/2.megaBOLT/Normal/Normal.csv";
			$cmd.= " --list $outdir/$sampid/2.megaBOLT/$dtype.csv";
			$cmd.= " --ref $configs{reference}" if($configs{reference});
			$cmd.= " --knownSites $configs{Cosmic}" if($configs{Cosmic});
			$cmd.= " --knownSites $configs{GATKdbsnp}" if($configs{GATKdbsnp});
			$cmd.= " --knownSites $configs{'1kg_phase1_snp'}" if($configs{'1kg_phase1_snp'});
			$cmd.= " --knownSites $configs{Mills_indel}" if($configs{Mills_indel});
			$cmd.= " --vcf $configs{GATKdbsnp}" if($configs{GATKdbsnp});
			if($configs{tumor_only} =~ /true/i && $dtype =~ /Tumor/i){
				$cmd.= " && \\\n";
				my $mydir = "$outdir/$sampid/5.somatic/tumor_only";
				system("mkdir -p $mydir") unless(-e $mydir);
				$cmd.= "MegaBOLT --type mutect2 --mutect2-input $bqsrbam{$sampid}{$dtype} --tumor $sampid-Tumor";
				$cmd.= " --ref $configs{reference}" if($configs{reference});
				$cmd.= " --vcf $configs{GATKdbsnp}" if($configs{GATKdbsnp});
				$cmd.= " --knownSites $configs{Cosmic}" if($configs{Cosmic});
				$cmd.= " --knownSites $configs{GATKdbsnp}" if($configs{GATKdbsnp});
				$cmd.= " --knownSites $configs{'1kg_phase1_snp'}" if($configs{'1kg_phase1_snp'});
				$cmd.= " --knownSites $configs{Mills_indel}" if($configs{Mills_indel});		
				$cmd.= " --outputdir $mydir && \\\n";
				$cmd.= "cp $mydir/output/output.mutect2.vcf.gz $mydir/$sampid.somatic.vcf.gz && \\\n";
				$cmd.= "cp $mydir/output/output.mutect2.vcf.gz.tbi $mydir/$sampid.somatic.vcf.gz.tbi && \\\n";
				$cmd.= "rm -rf $mydir/output";
			}
			if($configs{metaphlan} =~ /true/i || $configs{kraken} =~ /true/i || $configs{hlaminer} =~ /true/i){
				$cmd.= " && \\\n";
				$cmd.= "$configs{samtools} fastq -@ 10 -N -f 4 -o $unmapfqs{$sampid}{$dtype} $bqsrbam{$sampid}{$dtype}";
			}
			if($configs{hlaminer} =~ /true/i){
				$cmd.= " && \\\n";
				$cmd.= "$configs{samtools} view -h -@ 10 -L $configs{TASRbed} $bqsrbam{$sampid}{$dtype} |";
				$cmd.= " $configs{samtools} fastq -@ 10 -N -o $outdir/$sampid/2.megaBOLT/$sampid-$dtype/HLA.fq.gz";
			}
			generateShell($shell, $cmd);
			$scripts{MegaBOLT}{$sampid}{$dtype} = $shell;

			# germline mutation annotation
			$shell = "$outdir/$sampid/0.shell/a02.megaBOLT/a02.annosnp_$dtype.sh";
			$cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
			$cmd.= "rm $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz.tbi\n";
			$cmd.= "gzip -d $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.hc.vcf.gz\n";
			$cmd.= "$configs{GATK4} VariantFiltration -R $configs{reference} $configs{SNPfilterParameter}";
			$cmd.= " $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.hc.vcf";
			$cmd.= " --verbosity ERROR -O $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf && \\\n";
			$cmd.= "rm $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.hc.vcf && \\\n";
			$cmd.= "awk '\$1~/^#/ || \$7==\"PASS\"' $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf >";
			$cmd.= " $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf && \\\n";
			$cmd.= "$configs{bgzip} -c $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf >";
			$cmd.= " $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz && \\\n";
			$cmd.= "$configs{tabix} -p vcf $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz && \\\n";
			$cmd.= "$configs{bgzip} -@ 2 -f $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf && \\\n";
			$cmd.= "$configs{tabix} -p vcf $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf.gz && \\\n";
			$cmd.= "rm $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf && \\\n";
			$cmd.= "perl $configs{ANNOVAR} $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz";
			$cmd.= " $configs{ANNOVAR_refdir} -out $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.germline.variants.annot $configs{annovar_par} && \\\n";
			$cmd.= "#less $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.germline.variants.annot.homo_multianno.txt | cut -f1,2,4,5,6 |";
			$cmd.= " grep 'rs' > $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.final.variants.txt && \\\n";
			$cmd.= "ls $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.germline.variants.annot* | while read file;";
			$cmd.= "do $configs{bgzip} -@ 2 -f \$file; done";
			generateShell($shell,$cmd);
			$scripts{gatksnp_ann}{$sampid}{$dtype}=$shell;

			##snp chemicalDrug
			$shell = "$outdir/$sampid/0.shell/a16.medicine/a16.chemicalDrug_$dtype.sh";
			$cmd = "perl $configs{chemicaldrug_tool} $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.final.variants.txt >";
			$cmd.= " $outdir/$sampid/9.medicine/chemicalDrug/chemicalDrug.out";
			generateShell($shell,$cmd);
			$scripts{chemicalDurg_ann}{$sampid}{$dtype}=$shell;
			# }

			if($configs{tumor_only} =~ /true/i && $dtype =~ /Tumor/i){
				# tumor only result filter
				my $mydir = "$outdir/$sampid/5.somatic/tumor_only";
				system("mkdir -p $mydir/contamination/tmp") unless (-e "$mydir/contamination/tmp");
				system("mkdir -p $mydir/artifact/tmp") unless (-e "$mydir/artifact/tmp");
				my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/tumor_only.bqsr.contamination.sh";
				my $tbam = $bqsrbam{$sampid}{$dtype};
				my $cmd = "$configs{GATK4} --java-options \"-Xmx6g -XX:ParallelGCThreads=1";
				$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\"";
				$cmd.= " GetPileupSummaries -I $tbam -V $configs{GetPileupSummaries}";
				$cmd.= " -L $configs{reference_bed}";
				$cmd.= " -O $mydir/contamination/$sampid-Tumor.pileups.table && \\\n";
				$cmd.= "$configs{GATK4} --java-options \"-Xmx6g -XX:ParallelGCThreads=1";
				$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\"";
				$cmd.= " CalculateContamination -I $mydir/contamination/$sampid-Tumor.pileups.table";
				$cmd.= " -O $mydir/contamination/$sampid.contamination.table";
				generateShell($shell,$cmd);
				$scripts{tumor_only}{$sampid}{contamination} = $shell;
				
				## Artifact calculation
				# $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/tumor_only.bqsr.artifact.sh";
				# $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
				# $cmd.= "$configs{GATK4} --java-options \"-Xmx10g -XX:ParallelGCThreads=1";
				# $cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/artifact/tmp\"";
				# $cmd.= " CollectSequencingArtifactMetrics -R $configs{reference} -I $tbam --FILE_EXTENSION \".txt\"";
				# $cmd.= " -O $mydir/artifact/$sampid.artifact";
				# generateShell($shell,$cmd);
				# $scripts{tumor_only}{$sampid}{artifact} = $shell;

				## Filter somatic mutation Calls
				$shell = "$outdir/$sampid/0.shell/a02.megaBOLT/tumor_only.bqsr.filter.sh";
				$cmd = "mkdir -p $mydir/tmp && \\\n";
				$cmd.= "$configs{GATK4} --java-options \"-Xmx25g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
				$cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterMutectCalls";
				$cmd.= " -V $mydir/$sampid.somatic.vcf.gz -contamination-table $mydir/contamination/$sampid.contamination.table";
				$cmd.= " -O $mydir/$sampid.somatic.filter1.vcf.gz && \\\n";
				# $cmd.= "$configs{GATK4} --java-options \"-Xmx25g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
				# $cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterByOrientationBias";
				# $cmd.= " --artifact-modes \"G/T\" --artifact-modes \"C/T\" -V $mydir/$sampid.somatic.filter1.vcf.gz";
				# $cmd.= " -P $mydir/artifact/$sampid.artifact.pre_adapter_detail_metrics.txt";
				# $cmd.= " -O $mydir/$sampid.somatic.filter2.vcf.gz && \\\n";
				# $cmd.= "gzip -dc $mydir/$sampid.somatic.filter2.vcf.gz |";
				$cmd.= "gzip -dc $mydir/$sampid.somatic.filter1.vcf.gz |";
				$cmd.= " perl -ne \'if(/^#/){print \$_;}else{my \$FILTER=(split /\\t/)[6];if(\$FILTER=~\"PASS\"){print \$_;}}\' |";
				$cmd.= " gzip -c > $mydir/$sampid.somatic.PASS.vcf.gz && \\\n";
				$cmd.= "perl $configs{ANNOVAR} $mydir/$sampid.somatic.PASS.vcf.gz";
				$cmd.= " $configs{ANNOVAR_refdir} -out $mydir/$sampid.somatic.annot $configs{annovar_par} && \\\n";
				$cmd.= "ls $mydir/$sampid.somatic.annot* | while read file;";
				$cmd.= "do $configs{bgzip} -@ 2 -f \$file; done && \\\n";
				$cmd.= "rm -rf $mydir/artifact $mydir/contamination $mydir/$sampid.somatic.vcf.gz* $mydir/$sampid.somatic.filter1* $mydir/$sampid.somatic.filter2*";
				generateShell($shell,$cmd);
				$scripts{tumor_only}{$sampid}{filter} = $shell;
			}
		}
	}
}

sub megaboltqc{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my ($tbam, $nbam);
		for my $dtype(keys %{$hash{$sampid}}){
			$tbam = $hash{$sampid}{$dtype} if($dtype =~ /Tumor/i);
			$nbam = $hash{$sampid}{$dtype} if($dtype =~ /Normal/i);
			my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/$dtype.coverage.sh";
			my $mydir = "$outdir/$sampid/2.megaBOLT/coverage";
			system("mkdir -p $mydir/$dtype") unless (-e "$mydir/$dtype");
			my $cmd = "$configs{samtools} bedcov $configs{targetRegion} $hash{$sampid}{$dtype} > $mydir/$dtype/bedcov.txt && \\\n";
			$cmd.= "$configs{samtools} coverage $hash{$sampid}{$dtype} > $mydir/$dtype/coverage.txt && \\\n";
			$cmd.= "$configs{samtools} depth -a -o $mydir/$dtype/depth.txt $hash{$sampid}{$dtype} && \\\n";
			$cmd.= "$configs{bgzip} -@ 2 $mydir/$dtype/depth.txt";
			generateShell($shell,$cmd);
			$scripts{coverage}{$sampid}{$dtype} = $shell;          
		}
		unless($configs{tumor_only} =~ /true/i || ! -e $configs{qcvcf}){
			my $mydir = "$outdir/$sampid/2.megaBOLT/QC";
			system("mkdir -p $mydir/tmp") unless (-e "$mydir/tmp");
			my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/QC.sh";
			my $cmd = "$configs{python} $configs{qctool} --bam1 $nbam --bam2 $tbam";
			$cmd.= " --output $mydir/bam_matcher.report.txt --vcf $configs{qcvcf}";
			$cmd.= " --reference $configs{reference} --scratch-dir $mydir/tmp $configs{qcpar} && \\\n";
			$cmd.= "rm -rf $mydir/tmp";
			generateShell($shell,$cmd);
			$scripts{megaboltqc}{$sampid} = $shell;
		}
	}
}

sub metaphlan{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			my $shell = "$outdir/$sampid/0.shell/a24.metaphlan/a24.metaphlan_$dtype.sh";
			my $fastq = $hash{$sampid}{$dtype};
			my $cmd = "export PATH=\"/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/02.apps/miniconda3/envs/metaphlan/bin:\$PATH\" && \\\n";
			$cmd.= "mkdir -p $outdir/$sampid/11.microbiome/$dtype && \\\n";
			$cmd.= "$configs{metaphlan3} $fastq --nproc $configs{metaphlan_nproc} --input_type fastq";
			$cmd.= " --bowtie2out $outdir/$sampid/11.microbiome/$dtype/${sampid}_${dtype}_metagenome.bowtie2.bz2";
			$cmd.= " -o $outdir/$sampid/11.microbiome/$dtype/${sampid}_${dtype}_metaphlan_wavg_l.txt";
			$cmd.= " --index $configs{metaphlan_index} --bowtie2db $configs{bowtie2db}";
			$cmd.= " --bowtie2_exe $configs{bowtie2} --bowtie2_build $configs{bowtie2_build}";
			$cmd.= " --sample_id ${sampid}_${dtype} --add_viruses --stat wavg_l && \\\n";
			$cmd.= "$configs{metaphlan3} $outdir/$sampid/11.microbiome/$dtype/${sampid}_${dtype}_metagenome.bowtie2.bz2 --input_type bowtie2out";
			$cmd.= " -o $outdir/$sampid/11.microbiome/$dtype/${sampid}_${dtype}_metaphlan_wavg_g.txt";
			$cmd.= " --index $configs{metaphlan_index} --bowtie2db $configs{bowtie2db}";
			$cmd.= " --sample_id ${sampid}_${dtype} --add_viruses --stat wavg_g";
			generateShell($shell,$cmd);
			$scripts{metaphlan}{$sampid}{$dtype} = $shell;
		}
	}
}

sub kraken2{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			my $shell = "$outdir/$sampid/0.shell/a25.kraken2/a25.kraken2_$dtype.sh";
			my $fastq = $hash{$sampid}{$dtype};
			my $cmd = "$configs{kraken2} --db $configs{Kraken2DB} --threads $configs{Kraken2_threads} $configs{Kraken2par}";
			$cmd.= " --output $outdir/$sampid/11.microbiome/$dtype/${sampid}_${dtype}_kraken2.out";
			$cmd.= " --report $outdir/$sampid/11.microbiome/$dtype/${sampid}_${dtype}_kraken2.report $fastq";
			generateShell($shell, $cmd);
			$scripts{kraken2}{$sampid}{$dtype} = $shell;
		}
	}
}

sub Manta{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/3.SV/Manta";
		system("mkdir -p $mydir") unless (-e "$mydir");
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		my $shell = "$outdir/$sampid/0.shell/a20.Manta/Manta.sh";
		my $cmd = "$configs{python} $configs{Mantatool} --tumorBam $tbam --normalBam $nbam";
		$cmd.= " --referenceFasta $configs{reference} --runDir $mydir && \\\n";
		$cmd.= "$configs{python} $mydir/runWorkflow.py -m local -j 8";
		generateShell($shell,$cmd);
		$scripts{Manta}{$sampid} = $shell;
	}
}

sub svaba{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/svaba";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		for my $mychr(@chrs){
			my $shell = "$outdir/$sampid/0.shell/a12.svaba/svaba.$mychr.sh";
			my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
			$cmd.= "$configs{svabatool} run -G $configs{reference} -t $tbam -n $nbam -k $configs{reference_block_bed_path}/$mychr.bed";
			$cmd.= " -D $configs{GATKdbsnp}" if($configs{GATKdbsnp});
			$cmd.= " -a $mydir/$mychr";
			generateShell($shell,$cmd);
			$scripts{svaba}{$sampid}{$mychr} = $shell;
		}
		my $shell = "$outdir/$sampid/0.shell/a12.svaba/merge.svaba.sh";
		my $cmd = "cat $mydir/chr1.svaba.somatic.indel.vcf | grep '^##' | grep -v '^##source'> $mydir/svaba.head && \\\n";
		$cmd.= "cat $mydir/chr*.svaba.somatic.indel.vcf | grep '^##source' >> $mydir/svaba.head && \\\n";
		$cmd.= "cat $mydir/chr1.svaba.somatic.indel.vcf | grep '^#CHROM' >> $mydir/svaba.head && \\\n";
		$cmd.= "cat $mydir/chr*.svaba.somatic.indel.vcf | grep -v '#' | grep PASS | cat $mydir/svaba.head - > $mydir/$sampid.svaba.vcf && \\\n";
		$cmd.= "gzip -f $mydir/$sampid.svaba.vcf && \\\n";
		$cmd.= "rm $mydir/chr*";
		generateShell($shell,$cmd);
		$scripts{svaba}{$sampid}{merge} = $shell;
		$indelvcf{$sampid}{svaba}="$mydir/$sampid.svaba.vcf.gz";
	}
}

sub mutectSNV{
	my (%hash) = @_;    # BQSR bam, whole Tumor/Normal bam
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/mutect";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		for my $mychr(@chrs){
			my $shell = "$outdir/$sampid/0.shell/a05.mutect/mutect.bqsr.$mychr.sh";
			my $thistmp = "$mydir/$mychr.tmp";
			system("mkdir -p $thistmp") unless (-e "$thistmp");
			my $cmd = "rm -f $mydir/$sampid.mutect.$mychr.txt.gz\n";
			$cmd.= "$configs{java17} -Xmx10g -Djava.io.tmpdir=$thistmp -XX:-UseGCOverheadLimit";
			$cmd.= " -jar $configs{mutectjar} -T MuTect -R $configs{reference}";
			$cmd.= " --input_file:normal $nbam --input_file:tumor $tbam";
			$cmd.= " --normal_sample_name $sampid-Normal --tumor_sample_name $sampid-Tumor";
			$cmd.= " --vcf $mydir/$sampid.mutect.$mychr.vcf.gz --out $mydir/$sampid.mutect.$mychr.txt";
			$cmd.= " --dbsnp $configs{dbsnp}" if($configs{dbsnp});
			$cmd.= " --cosmic $configs{Cosmic}" if($configs{Cosmic});
			# $cmd.= " -L $outdir/$sampid/TR/Splited_by_chr/$mychr.bed";
			$cmd.= " -L $configs{reference_block_bed_path}/$mychr.bed";
			$cmd.= " --enable_extended_output --downsampling_type NONE && \\\n";
			$cmd.= "gzip $mydir/$sampid.mutect.$mychr.txt";
			generateShell($shell,$cmd);
			$scripts{mutect}{step1}{$sampid}{$mychr} = $shell;
		}
		my $shell = "$outdir/$sampid/0.shell/a05.mutect/mergevcf.mutect.sh";
		my $cmd = "gzip -dc $mydir/$sampid.mutect.chr1.txt.gz | head -2 > $mydir/head.txt && \\\n";
		$cmd.= "gzip -dc $mydir/$sampid.mutect.*.txt.gz | grep \"$sampid\" |";
		$cmd.= " cat $mydir/head.txt - | gzip -f > $mydir/MuTect.txt.gz && \\\n";
		$cmd.= "gzip -dc $mydir/$sampid.mutect.chr1.vcf.gz | grep '^#' > $mydir/head2.txt && \\\n";
		$cmd.= "gzip -dc $mydir/$sampid.mutect.*.vcf.gz | grep -v '^#' |";
		$cmd.= " cat $mydir/head2.txt - | gzip -f > $mydir/$sampid.MuTect.vcf.gz && \\\n";
		# $cmd.= "perl $configs{mutectFilter} $mydir/MuTect.txt.gz $mydir/$sampid.MuTect.vcf.gz 0.02 0.05";
		# $cmd.= " $mydir/$sampid.mutect.merge.txt.gz $mydir/$sampid.mutect.merge.vcf.gz";
		$cmd.= "rm -rf $mydir/*.tmp && \\\n";
		$cmd.= "rm $mydir/head* && \\\n";
		$cmd.= "rm $mydir/$sampid.mutect.*.txt.gz && \\\n";
		$cmd.= "rm $mydir/$sampid.mutect.*.vcf.gz && \\\n";
		$cmd.= "rm $mydir/$sampid.mutect.*.vcf.gz.idx";
		generateShell($shell,$cmd);
		$scripts{mutect}{step2}{$sampid} = $shell;
		# $snvvcf{$sampid}{mutect} = "$mydir/$sampid.mutect.merge.vcf.gz";
		$snvvcf{$sampid}{mutect} = "$mydir/$sampid.MuTect.vcf.gz";
		# $snvvcf{$sampid}{mutectraw} = "$mydir/$sampid.MuTect.vcf.gz";
	}
}

sub mutect2{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/mutect2";
		system("mkdir -p $mydir/contamination/tmp") unless (-e "$mydir/contamination/tmp");
		system("mkdir -p $mydir/artifact/tmp") unless (-e "$mydir/artifact/tmp");
		system("mkdir -p $mydir/PoN/tmp") unless (-e "$mydir/PoN/tmp");
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		## Contamination Estimate
		my $shell = "$outdir/$sampid/0.shell/a06.mutect2/mutect2.bqsr.contamination.sh";
		my $cmd = "$configs{GATK4} --java-options \"-Xmx6g -XX:ParallelGCThreads=1";
		$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\"";
		$cmd.= " GetPileupSummaries -I $tbam -V $configs{GetPileupSummaries}";
		$cmd.= " -L $configs{reference_bed}";
		$cmd.= " -O $mydir/contamination/$sampid-Tumor.pileups.table && \\\n";
		$cmd.= "$configs{GATK4} --java-options \"-Xmx6g -XX:ParallelGCThreads=1";
		$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\"";
		$cmd.= " GetPileupSummaries -I $nbam -V $configs{GetPileupSummaries}";
		$cmd.= " -L $configs{reference_bed}";
		$cmd.= " -O $mydir/contamination/$sampid-Normal.pileups.table && \\\n";
		$cmd.= "$configs{GATK4} --java-options \"-Xmx4g -XX:ParallelGCThreads=1";
		$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\"";
		$cmd.= " CalculateContamination -I $mydir/contamination/$sampid-Tumor.pileups.table";
		$cmd.= " -matched $mydir/contamination/$sampid-Normal.pileups.table";
		$cmd.= " -O $mydir/contamination/$sampid.contamination.table";
		generateShell($shell,$cmd);
		$scripts{mutect2}{$sampid}{contamination} = $shell;
		## Artifact calculation
		$shell = "$outdir/$sampid/0.shell/a06.mutect2/mutect2.bqsr.artifact.sh";
		$cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
		$cmd.= "$configs{GATK4} --java-options \"-Xmx4g -XX:ParallelGCThreads=1";
		$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/artifact/tmp\"";
		$cmd.= " CollectSequencingArtifactMetrics -R $configs{reference} -I $tbam --FILE_EXTENSION \".txt\"";
		$cmd.= " -O $mydir/artifact/$sampid.artifact";
		generateShell($shell,$cmd);
		$scripts{mutect2}{$sampid}{artifact} = $shell;
		##normal Panel
		$shell = "$outdir/$sampid/0.shell/a06.mutect2/CreateSomaticPanelOfNormals.sh";
		$cmd = "mkdir -p $mydir/PoN/tmp\n";
		$cmd.= "$configs{GATK4} --java-options \"-Xmx16g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
		$cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/PoN/tmp\" Mutect2 -R $configs{reference}";
		$cmd.= " -I $nbam -tumor $sampid-Normal --germline-resource $configs{gnomAD}";
		$cmd.= " -O $mydir/PoN/$sampid-Normal.vcf.gz -L $configs{targetRegion}";
		generateShell($shell,$cmd) unless($configs{panel_of_normal});
		$scripts{mutect2}{$sampid}{PoN} = $shell;
		## Mutect2
		my @vcfs;
		# my @splitbeds=glob "$outdir/$sampid/TR/Splited_by_chr/chr*.bed";
		# for my $mybed(@splitbeds){
			# my $mychr=basename($mybed);
			# $mychr=~s/\.bed//gi;
		for my $mychr(@chrs){
			my $shell = "$outdir/$sampid/0.shell/a06.mutect2/mutect2.bqsr.call.$mychr.sh";
			my $thistmp = "$mydir/$mychr.tmp";
			system("mkdir -p $thistmp") unless (-e "$thistmp");
			$cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
			$cmd.= "$configs{GATK4} --java-options \"-Xmx8g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
			$cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$thistmp\" Mutect2";
			$cmd.= " -R $configs{reference} -I $tbam -tumor $sampid-Tumor -I $nbam -normal $sampid-Normal";
			$cmd.= " --germline-resource $configs{gnomAD} --af-of-alleles-not-in-resource 0.00003125";
			if($configs{panel_of_normal}){
				$cmd.= " --panel-of-normals $configs{panel_of_normal}";
			}
			else{
				$cmd.= " --panel-of-normals $mydir/PoN/$sampid-Normal.vcf.gz";
			}
			$cmd.= " -L $configs{reference_block_bed_path}/$mychr.bed";
			$cmd.= " -O $mydir/$sampid.mutect2.$mychr.vcf.gz";
			generateShell($shell,$cmd);
			$scripts{mutect2}{step1}{$sampid}{$mychr} = $shell;
			push @vcfs, "$mydir/$sampid.mutect2.$mychr.vcf.gz";
		}
		$shell = "$outdir/$sampid/0.shell/a06.mutect2/mergevcf.mutect2.sh";
		$cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
		$cmd.= "$configs{vcfconcat} @vcfs | $configs{vcfsort} -t $mydir > $mydir/$sampid.mutect2.vcf && \\\n";
		$cmd.= "$configs{bgzip} $mydir/$sampid.mutect2.vcf && \\\n";
		$cmd.= "$configs{tabix} $mydir/$sampid.mutect2.vcf.gz";
		generateShell($shell,$cmd);
		$scripts{mutect2}{step2}{$sampid} = $shell;
		## Filter Mutect Calls
		$shell = "$outdir/$sampid/0.shell/a06.mutect2/mutect2.bqsr.filter.sh";
		$cmd = "mkdir -p $mydir/tmp && \\\n";
		$cmd.= "$configs{GATK4} --java-options \"-Xmx25g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
		$cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterMutectCalls";
		$cmd.= " -V $mydir/$sampid.mutect2.vcf.gz -contamination-table $mydir/contamination/$sampid.contamination.table";
		$cmd.= " -O $mydir/$sampid.mutect2.filter1.vcf.gz && \\\n";
		$cmd.= "$configs{GATK4} --java-options \"-Xmx25g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
		$cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterByOrientationBias";
		$cmd.= " --artifact-modes \"G/T\" --artifact-modes \"C/T\" -V $mydir/$sampid.mutect2.filter1.vcf.gz";
		$cmd.= " -P $mydir/artifact/$sampid.artifact.pre_adapter_detail_metrics.txt";
		$cmd.= " -O $mydir/$sampid.mutect2.filter2.vcf.gz && \\\n";
		$cmd.= "rm -rf $mydir/chr*.tmp $mydir/$sampid.mutect2.chr*";
		generateShell($shell,$cmd);
		$scripts{mutect2}{$sampid}{filter} = $shell;
		$snvvcf{$sampid}{mutect2} = "$mydir/$sampid.mutect2.filter2.vcf.gz";
		$indelvcf{$sampid}{mutect2} = "$mydir/$sampid.mutect2.filter2.vcf.gz";
	}    
}

=removed
sub msi{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/msi";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		my $shell = "$outdir/$sampid/0.shell/a21.msi/mantis.sh";
		my $cmd = "$configs{python3} $configs{mantis_tool} -t $tbam -n $nbam -b $configs{msi_bed}";
		$cmd.= " $configs{mantispar} --genome $configs{mantis_reference} -o $mydir/mantis.msi.txt";
		generateShell($shell,$cmd);
		$scripts{mantis}{$sampid} = $shell;
		$shell= "$outdir/$sampid/0.shell/a21.msi/tumor.reheader.sh";
		$cmd = "$configs{samtools} view -H $tbam >$mydir/tumor.header.sam && \\\n";
		$cmd.= "less -S $mydir/tumor.header.sam | grep -v '\@CO' > $mydir/tumor.reheader.sam && \\\n";
		$cmd.= "$configs{samtools} reheader -P $mydir/tumor.reheader.sam $tbam >$mydir/tumor.bam && \\\n";
		$cmd.= "$configs{samtools} index -b $mydir/tumor.bam";
		generateShell($shell,$cmd);
		$scripts{tumor_reheader}{$sampid} = $shell;
		$shell= "$outdir/$sampid/0.shell/a21.msi/normal.reheader.sh";
		$cmd = "$configs{samtools} view -H $nbam >$mydir/normal.header.sam && \\\n";
		$cmd.= "less -S $mydir/normal.header.sam | grep -v '\@CO' > $mydir/normal.reheader.sam && \\\n";
		$cmd.= "$configs{samtools} reheader -P $mydir/normal.reheader.sam $nbam >$mydir/normal.bam && \\\n";
		$cmd.= "$configs{samtools} index -b $mydir/normal.bam";
		generateShell($shell,$cmd);
		$scripts{normal_reheader}{$sampid} = $shell;
		$shell= "$outdir/$sampid/0.shell/a21.msi/missensor.sh";
		$cmd = "$configs{msisensor_tool} msi -d $configs{msisor_bed} -t $mydir/tumor.bam -n $mydir/normal.bam";
		$cmd.= " -o $mydir/msisensor.txt -l 1 -q 3 && \\\n";
		$cmd.= "rm $mydir/tumor.bam $mydir/normal.bam";
		generateShell($shell,$cmd);
		$scripts{msisensor}{$sampid} = $shell;
	}
}
=cut

sub msi{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/msi";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		my $shell = "$outdir/$sampid/0.shell/a21.msi/mantis.sh";
		my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu && \\\n";
		$cmd.= "$configs{python3} $configs{mantis_tool} -t $tbam -n $nbam -b $configs{msi_bed} $configs{mantispar}";
		$cmd.= " --genome $configs{mantis_reference} --threads 10 -o $mydir/$sampid.mantis.msi.txt";
		generateShell($shell,$cmd);
		$scripts{mantis}{$sampid} = $shell;

		$shell= "$outdir/$sampid/0.shell/a21.msi/missensor.sh";
		$cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n"; # [Ming-20191129] source
		$cmd.= "$configs{msisensor_tool} msi -d $configs{msisor_bed} -t $tbam -n $nbam -o $mydir/$sampid.msisensor.txt";
		generateShell($shell,$cmd);
		$scripts{msisensor}{$sampid} = $shell;
	}
}

sub strelka{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/strelka";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		my $shell = "$outdir/$sampid/0.shell/a07.strelka/strelka.$sampid.sh";
		my $cmd = "rm -rf $mydir \n";
		$cmd.= "$configs{strelkapl} --normal=$nbam --tumor=$tbam --ref=$configs{reference}";
		$cmd.= " --config=$configs{strelkapar} --output-dir=$mydir";
		generateShell($shell, $cmd);
		$scripts{strelka}{step1}{$sampid} = $shell;
		for my $mychr(@chrs){
			$shell = "$outdir/$sampid/0.shell/a07.strelka/runstrelka.$mychr.sh";
			$cmd = "if [ -d \"$mydir/$mychr\" ];then\n";
			$cmd.= "\tcd $mydir/$mychr && make\n";
			$cmd.= "fi";
			generateShell($shell, $cmd);
			$scripts{strelka}{step2}{$sampid}{$mychr} = $shell;
		}
		$shell = "$outdir/$sampid/0.shell/a07.strelka/mergeStrelka.sh";
		$cmd = "less $mydir/chr1/results/passed.somatic.snvs.vcf | grep '^#' > $mydir/head.snv.txt && \\\n";
		$cmd.= "less $mydir/chr1/results/passed.somatic.indels.vcf | grep '#' > $mydir/head.indel.txt && \\\n";
		$cmd.= "cat $mydir/*/results/passed.somatic.snvs.vcf | grep -v '#' | grep 'PASS' | cat $mydir/head.snv.txt - |";
		$cmd.= " gzip -f > $mydir/$sampid.strelka.snv.vcf.gz && \\\n";
		$cmd.= "cat $mydir/*/results/passed.somatic.indels.vcf | grep -v '#' | grep 'PASS' | cat $mydir/head.indel.txt - |";
		$cmd.= " gzip -f > $mydir/$sampid.strelka.indel.vcf.gz && \\\n";
		$cmd.= "rm -r $mydir/chr* $mydir/head.snv.txt $mydir/head.indel.txt";
		generateShell($shell,$cmd);
		$scripts{strelka}{step3}{$sampid}=$shell;        
		$snvvcf{$sampid}{strelka} = "$mydir/$sampid.strelka.snv.vcf.gz";
		$indelvcf{$sampid}{strelka} = "$mydir/$sampid.strelka.indel.vcf.gz";
	}
}

sub strelka2{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/strelka2";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		my $shell = "$outdir/$sampid/0.shell/a08.strelka2/$sampid.strelks2.sh";
		my $cmd = "rm -rf $mydir/*\n";
		$cmd.= "$configs{python} $configs{strelka2py} --normalBam $nbam --tumorBam $tbam";
		# $cmd.= " --referenceFasta $configs{reference} --runDir $mydir --callRegions $configs{bedgzip} && \\\n";
		$cmd.= " --referenceFasta $configs{reference} --runDir $mydir && \\\n";
		$cmd.= "$configs{python} $mydir/runWorkflow.py $configs{runstrelka2flow} && \\\n";
		$cmd.= "ln -s $mydir/results/variants/somatic.snvs.vcf.gz $mydir/$sampid.strelks2.snv.vcf.gz && \\\n";
		$cmd.= "ln -s $mydir/results/variants/somatic.indels.vcf.gz $mydir/$sampid.strelks2.indel.vcf.gz";
		generateShell($shell,$cmd);
		$scripts{strelka2}{$sampid} = $shell;
		$snvvcf{$sampid}{strelka2} = "$mydir/$sampid.strelks2.snv.vcf.gz";
		$indelvcf{$sampid}{strelka2} = "$mydir/$sampid.strelks2.indel.vcf.gz";
	}
}

sub varscan2{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/varscan2";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		my @snvvcfs; my @indelvcfs;
		for my $mychr(@chrs){
			my $shell = "$outdir/$sampid/0.shell/a22.varscan2/varscan2.bqsr.$mychr.sh";
			my $cmd = "mkfifo $mydir/$mychr.fifo\n";
			$cmd.= "$configs{samtools} mpileup -q 20 -Q 20 -B";
			# $cmd.= " -l $outdir/$sampid/TR/Splited_by_chr/$mychr.bed";
			$cmd.= " -l $configs{reference_block_bed_path}/$mychr.bed";
			$cmd.= " -f $configs{reference} $nbam $tbam > $mydir/$mychr.fifo &";
			$cmd.= " $configs{java} -Xmx3g -jar $configs{varscan2jar}";
			$cmd.= " somatic $mydir/$mychr.fifo --output-snp $mydir/varscan2.$mychr.snv";
			$cmd.= " --output-indel $mydir/varscan2.$mychr.indel --min-var-freq 0.05 --mpileup 1";
			$cmd.= " --min-freq-for-hom 0.75 --somatic-p-value 0.05 --p-value 0.1 --min-coverage 3 --strand-filter 1 && \\\n";
			$cmd.= "rm -f $mydir/$mychr.fifo";
			generateShell($shell,$cmd);
			# $scripts{varscan2}{$sampid} = $shell;
			$scripts{varscan2}{step1}{$sampid}{$mychr} = $shell;
			}
		my $shell = "$outdir/$sampid/0.shell/a22.varscan2/varscan2.mergevcf.sh";
		my $cmd = "grep -v 'chrom' $mydir/varscan2.chr*.snv > $mydir/$sampid.varscan2.merge.content.snv && \\\n";
		$cmd.= "grep -v 'chrom' $mydir/varscan2.chr*.indel > $mydir/$sampid.varscan2.merge.content.indel && \\\n";
		$cmd.= "head -1 $mydir/varscan2.chr1.snv >$mydir/$sampid.varscan2.head.txt && \\\n";
		$cmd.= "cat $mydir/$sampid.varscan2.head.txt $mydir/$sampid.varscan2.merge.content.snv >";
		$cmd.= "$mydir/$sampid.varscan2.merge.raw.snv && \\\n";
		$cmd.= "cat $mydir/$sampid.varscan2.head.txt $mydir/$sampid.varscan2.merge.content.indel >";
		$cmd.= "$mydir/$sampid.varscan2.merge.raw.indel && \\\n";
		$cmd.= "$configs{java} -Xmx3g -jar $configs{varscan2jar} processSomatic";
		$cmd.= " $mydir/$sampid.varscan2.merge.raw.snv --min-tumor-freq 0.05 --max-normal-freq 0.05 --p-value 0.05 && \\\n";
		$cmd.= "$configs{java} -Xmx3g -jar $configs{varscan2jar} processSomatic";
		$cmd.= " $mydir/$sampid.varscan2.merge.raw.indel --min-tumor-freq 0.05 --max-normal-freq 0.05 --p-value 0.05 && \\\n";
		$cmd.= "$configs{java} -Xmx3g -jar $configs{varscan2jar} somaticFilter";
		$cmd.= " $mydir/$sampid.varscan2.merge.raw.snv.Somatic --min-var-freq 0.05 --min-strands2 1 --p-value 0.05";
		$cmd.= " --indel-file $mydir/$sampid.varscan2.merge.raw.indel --output-file $mydir/$sampid.varscan2.merge.raw2.snv && \\\n";
		$cmd.= "$configs{java} -Xmx3g -jar $configs{varscan2jar} somaticFilter";
		$cmd.= " $mydir/$sampid.varscan2.merge.raw.indel.Somatic --min-var-freq 0.05 --min-strands2 1 --p-value 0.05";
		$cmd.= " --output-file $mydir/$sampid.varscan2.merge.raw2.indel && \\\n";
		$cmd.= "less $mydir/$sampid.varscan2.merge.raw2.snv |awk '\$18>0 && \$19>0' >$mydir/$sampid.varscan2.merge.snv.xls && \\\n";
		$cmd.= "less $mydir/$sampid.varscan2.merge.raw2.indel|awk '\$18>0 && \$19>0' >$mydir/$sampid.varscan2.merge.indel.xls && \\\n";
		$cmd.= "rm -rf $mydir/varscan2.chr*.snv $mydir/varscan2.chr*.indel && \\\n";
		$cmd.= "rm -rf $mydir/$sampid.varscan2.merge.raw2.snv $mydir/$sampid.varscan2.merge.raw2.indel && \\\n";
		$cmd.= "rm -rf $mydir/$sampid.varscan2.head.txt $mydir/$sampid.varscan2.merge.raw.snv $mydir/$sampid.varscan2.merge.raw.indel";
		#$cmd.= "$configs{bgzip} $mydir/$sampid.varscan2.snv.vcf && \\\n";
		#$cmd.= "$configs{tabix} $mydir/$sampid.varscan2.snv.vcf.gz && \\\n";
		generateShell($shell,$cmd);
		$scripts{varscan2}{step2}{$sampid} = $shell;
		$snvvcf{$sampid}{varscan2} = "$mydir/$sampid.varscan2.merge.snv.xls";
		$indelvcf{$sampid}{varscan2} = "$mydir/$sampid.varscan2.merge.indel.xls";
	}   
}

sub MuSE{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
	my $mydir = "$outdir/$sampid/5.somatic/MuSE";
	my $tbam = $hash{$sampid}{Tumor};
	my $nbam = $hash{$sampid}{Normal};
	# my @vcfs;
	for my $mychr(@chrs){
		my $shell = "$outdir/$sampid/0.shell/a23.MuSE/MuSE.bqsr.$mychr.sh";
		my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
		# $cmd.= "/ldfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\ \n";
		$cmd.= "$configs{muse} call -f $configs{reference}";
		$cmd.= " -l $configs{reference_block_bed_path}/$mychr.bed";
		$cmd.= " -O $mydir/$mychr $tbam $nbam";
		# $cmd.= "$configs{muse} sump -I $mydir/$sampid.$mychr.MuSE.txt -O $mydir/$sampid.$mychr.MuSE.vcf -D $configs{GATKdbsnp} -G";
		generateShell($shell,$cmd);
		$scripts{MuSE}{step1}{$sampid}{$mychr} = $shell;
		# push @vcfs,"$mydir/$mychr.MuSE.vcf";
	}
	my $shell = "$outdir/$sampid/0.shell/a23.MuSE/MuSE.mergevcf.sh";
	# my $cmd = "$configs{vcfconcat} @vcfs > $mydir/$sampid.MuSE.snv.raw.vcf && \\\n";
	# $cmd.= "less $mydir/$sampid.MuSE.snv.raw.vcf |grep -v 'Tier5' > $mydir/$sampid.MuSE.snv.vcf && \\\n";
	# $cmd.= "rm -rf $mydir/$sampid.MuSE.snv.raw.vcf && \\\n";
	# $cmd.= "$configs{bgzip} $mydir/$sampid.MuSE.snv.vcf && \\\n";
	# $cmd.= "$configs{tabix} $mydir/$sampid.MuSE.snv.vcf.gz";
	my $cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
	$cmd.= "less $mydir/chr1.MuSE.txt | grep '^#'> $mydir/txt.head.txt && \\\n";
	$cmd.= "cat $mydir/chr*.MuSE.txt | grep -v '^#' | cat $mydir/txt.head.txt - > $mydir/all.Muse.txt && \\\n";
	$cmd.= "$configs{muse} sump -I $mydir/all.Muse.txt -G -O $mydir/all.MuSE.vcf";
	$cmd.= " -D $configs{GATKdbsnp}" if($configs{GATKdbsnp});
	$cmd.= " && \\\n";
	$cmd.= "less $mydir/all.MuSE.vcf | grep '^#' > $mydir/head.txt && \\\n";
	$cmd.= "cat $mydir/all.MuSE.vcf | grep -v '^#' | grep 'PASS'| cat $mydir/head.txt - | gzip -f > $mydir/$sampid.MuSE.vcf.gz && \\\n";
	$cmd.= "rm $mydir/*MuSE.txt $mydir/head.txt $mydir/txt.head.txt ";
	generateShell($shell,$cmd);
	$scripts{MuSE}{step2}{$sampid} = $shell;
	$snvvcf{$sampid}{muse}="$mydir/$sampid.MuSE.vcf.gz";
  }
}

sub mergeSnvInDel{
	for my $sampid(sort keys %snvvcf){
		my $mydir = "$outdir/$sampid/5.somatic/result";
		my $mutect = $snvvcf{$sampid}{mutect} if($configs{mutect} =~ /true/i);
		my $muse = $snvvcf{$sampid}{muse} if($configs{MuSE} =~ /true/i);
		my $mutect2 = $snvvcf{$sampid}{mutect2} if($configs{mutect2} =~ /true/i);
		my $strelka = $snvvcf{$sampid}{strelka} if($configs{strelka} =~ /true/i);
		my $strelka_indel = $indelvcf{$sampid}{strelka} if($configs{strelka} =~ /true/i);
		my $strelka2 = $snvvcf{$sampid}{strelka2} if($configs{strelka2} =~ /true/i);       
		my $strelka2_indel = $indelvcf{$sampid}{strelka2} if($configs{strelka2} =~ /true/i);
		my $svaba = $indelvcf{$sampid}{svaba} if($configs{svaba} =~ /true/i);
		for my $mychr(@chrs){
			system("mkdir -p $mydir/tmp") unless (-e "$mydir/tmp");
			my $shell = "$outdir/$sampid/0.shell/a14.mergeMut/$sampid.merge_$mychr.sh";
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
		my $shell = "$outdir/$sampid/0.shell/a14.mergeMut/$sampid.concat.sh"; 
		my $cmd = "less $mydir/tmp/chr*.merged_snv.unsorted.vcf > $mydir/merged_snv.unsorted.vcf && \\\n";
		$cmd.= "less $mydir/tmp/chr*.uniq_snv.unsorted.vcf > $mydir/uniq_snv.unsorted.vcf && \\\n";
		$cmd.= "less $mydir/tmp/chr*.merged_indel.unsorted.vcf > $mydir/merged_indel.unsorted.vcf && \\\n";
		$cmd.= "less $mydir/tmp/chr*.uniq_indel.unsorted.vcf > $mydir/uniq_indel.unsorted.vcf && \\\n";
		$cmd.= "$configs{vcfsort} $mydir/merged_snv.unsorted.vcf | gzip -c > $mydir/$sampid.merged_snv.vcf.gz && \\\n";
		$cmd.= "$configs{vcfsort} $mydir/uniq_snv.unsorted.vcf | gzip -c > $mydir/$sampid.uniq_snv.vcf.gz && \\\n";
		$cmd.= "$configs{vcfsort} $mydir/merged_indel.unsorted.vcf | gzip -c > $mydir/$sampid.merged_indel.vcf.gz && \\\n";
		$cmd.= "$configs{vcfsort} $mydir/uniq_indel.unsorted.vcf | gzip -c > $mydir/$sampid.uniq_indel.vcf.gz && \\\n";
		$cmd.= "rm -rf $mydir/tmp $mydir/*.unsorted.vcf";
		generateShell($shell,$cmd);
		$scripts{mergeSnvInDel}{$sampid}{concat} = $shell;
		$shell = "$outdir/$sampid/0.shell/a14.mergeMut/phase.sh";
		$cmd = "perl $configs{ANNOVAR} $mydir/$sampid.merged_snv.vcf.gz $configs{ANNOVAR_refdir} -out $mydir/$sampid.merged_snv.annot $configs{annovar_par} && \\\n";
		$cmd.= "ls $mydir/$sampid.merged_snv.annot* | while read file; do $configs{bgzip} \$file; done && \\\n";
		$cmd.= "perl $configs{ANNOVAR} $mydir/$sampid.merged_indel.vcf.gz $configs{ANNOVAR_refdir} -out $mydir/$sampid.merged_indel.annot $configs{annovar_par} && \\\n";
		$cmd.= "ls $mydir/$sampid.merged_indel.annot* | while read file; do $configs{bgzip} \$file; done";
		$scripts{phase}{$sampid} = $shell;
		generateShell($shell,$cmd);
		$snvvcf{$sampid}{merge} = "$mydir/$sampid.merged_snv.vcf.gz";
		$indelvcf{$sampid}{merge} = "$mydir/$sampid.merged_indel.vcf.gz";
	}     
}

sub RNAbamfileterSNV{
	my ($vcf,$bam,$outfile) = @_;
	my $cmd = "perl $configs{checkBampl} -s $vcf -i $bam -o $outfile";
	return $cmd;
}

sub targetDrug{
	for my $sampid(sort keys %snvvcf){
		my $mutectvcf = $snvvcf{$sampid}{mutect};
		my $mutect2vcf = $snvvcf{$sampid}{mutect2};
		my $strelka2 = $indelvcf{$sampid}{strelka2};
		my $shell="$outdir/$sampid/0.shell/a16.medicine/Target_Drug_analysis.sh";
		#my $cmd = "perl $configs{target_bin}/SNV_onverlap.pl $snvvcf{$sampid}{merge}";
		# $cmd.= " $mutect2vcf $outdir/$sampid/9.medicine/targetDrug/Somatic.snv.vcf && \\\n";
		# $cmd.= "perl $configs{target_bin}/Sindel_onverlap.pl $indelvcf{$sampid}{merge} $mutect2vcf";
		# $cmd.= " $outdir/$sampid/9.medicine/targetDrug/Somatic.indel.vcf && \\\n";
		# $cmd.= "perl $configs{target_bin}/Sindel_filter.pl $outdir/$sampid/9.medicine/targetDrug/Somatic.indel.vcf";
		# $cmd.= " $outdir/$sampid/9.medicine/targetDrug/Somatic.indel.filter.vcf $sampid-Tumor 0.05 10 $configs{purity} && \\\n";
		# $cmd.= "perl $configs{target_bin}/SNV_Sindel_merge.pl $outdir/$sampid/9.medicine/targetDrug/Somatic.snv.vcf";
		# $cmd.= " $outdir/$sampid/9.medicine/targetDrug/Somatic.indel.filter.vcf";
		# $cmd.= " $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.vcf && \\\n";
		my $cmd ="less $outdir/$sampid/8.neoantigen/snv/Mutation/snv.annot.input | egrep 'splice|stop|frame|start_lost|missense_variant' |grep -v '^#'> $outdir/$sampid/8.neoantigen/snv/Mutation/snv.annot.input.filter && \\\n";
		$cmd.="less $outdir/$sampid/8.neoantigen/indel/Mutation/indel.annot.input | egrep 'splice|stop|frame|start_lost|missense_variant' |grep -v '^#'> $outdir/$sampid/8.neoantigen/indel/Mutation/indel.annot.input.filter && \\\n";
		$cmd.="cat $outdir/$sampid/8.neoantigen/snv/Mutation/snv.annot.input.filter $outdir/$sampid/8.neoantigen/indel/Mutation/indel.annot.input.filter >$outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.vcf && \\\n";
		$cmd.="sh $configs{vcf2maf} $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.vcf $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.maf $sampid-Normal $sampid-Tumor && \\\n";
		$cmd.="#add $sampid-Tumor to $configs{clinical} file \n";
		$cmd.="#run MafAnnotator.py at software-install node \n";
		$cmd.="#$configs{python} $configs{MafAnnotator} -i $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.maf -o $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.oncokb.maf -c $configs{clinical} \n";
		$cmd.="#perl $configs{tab2csv} $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.oncokb.maf > $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.oncokb.csv \n";
		$cmd.="gzip $outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.vcf";
		generateShell($shell,$cmd);    
		$scripts{targetDrug}{$sampid}=$shell;
		$snvvcf{$sampid}{target_maf}= "$outdir/$sampid/9.medicine/targetDrug/Somatic_SNV_indel.maf";
	}
}

sub FACETS{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
		my $mydir = "$outdir/$sampid/4.SCNV/FACETS";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};

		my $shell = "$outdir/$sampid/0.shell/a09.FACETS/FACETS.pileup.sh";
		my $cmd ="export LD_LIBRARY_PATH=/ldfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
		$cmd.= "if [ -e \"$mydir/pileup.gz\" ];then rm -rf $mydir/pileup.gz;fi && \\\n";
		$cmd.= "$configs{snppileup} $configs{FACETSpar} $configs{FACETSref} $mydir/pileup.gz $nbam $tbam";
		generateShell($shell,$cmd);
		$scripts{FACETS_pileup}{$sampid} = $shell;

		$shell = "$outdir/$sampid/0.shell/a09.FACETS/FACETS.calling.sh";
		$cmd = "export R_LIBS=/ldfssz1/ST_CANCER/POL/USER/lifuqiang/tools/library/R-3.3.2:";
		$cmd.= "/share/app/R-3.3.2/lib64/R/library && \\\n";
		$cmd.= "/share/app/R-3.5.2/bin/R $configs{FACETSTOOL} 9922 $mydir/pileup.gz";
		$cmd.= " $mydir/purity.txt $mydir/cnv.txt $mydir/cnv.pdf && \\\n";
		$cmd.= "gzip -f $mydir/cnv.txt";
		generateShell($shell,$cmd);
		$scripts{FACETS_call}{$sampid} = $shell;

		$shell = "$outdir/$sampid/0.shell/a09.FACETS/FACETS.anno.sh";
		$cmd = "zgrep -v \"^chrom\"  $mydir/cnv.txt.gz|";
		$cmd.= " perl -wlane 'print join \"\\t\", \"chr\$F[0]\", \$F[9], \$F[10], \"0\", \"0\", \"seg=\$F[1]\",";
		$cmd.= " \"num.mark=\$F[2]\", \"nhet=\$F[3]\", \"cnlr.median=\$F[4]\", \"mafR=\$F[5]\", \"segclust=\$F[6]\",";
		$cmd.= " \"cnlr.median.clust=\$F[7]\", \"mafR.clust=\$F[8]\", \"cf.em=\$F[11]\", \"tcn.em=\$F[12]\",";
		$cmd.= " \"lcn.em=\$F[13]\"' > $mydir/cnv.avinput && \\\n";
		$cmd.= "perl $configs{ANNOVAR} $mydir/cnv.avinput $configs{ANNOVAR_refdir}";
		$cmd.= " -buildver hg19 -out $mydir/cnv -remove";
		$cmd.= " -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,targetScanS,wgRna";
		$cmd.= " -operation g,g,g,r,r,r,r -nastring . && \\\n";
		$cmd.= "gzip -f $mydir/cnv.hg19_multianno.txt";
		generateShell($shell,$cmd);
		$scripts{FACETS_anno}{$sampid} = $shell;
	}
}

sub gatkcnv{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
		my $mydir = "$outdir/$sampid/4.SCNV/gatkcnv";
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		my $shell = "$outdir/$sampid/0.shell/a19.gatkcnv/gatkcnv.pileup.sh";
	
		`mkdir $mydir/GATKCNV` unless (-e "$mydir/GATKCNV");
		`mkdir $mydir/GATKCNV/config` unless (-e "$mydir/GATKCNV/config");
		`mkdir $mydir/GATKCNV/shell` unless (-e "$mydir/GATKCNV/shell");
		`mkdir $mydir/GATKCNV/result` unless (-e "$mydir/GATKCNV/result");

		open OUT,">$mydir/GATKCNV/config/bamlst" or die "$!";
		print OUT "$sampid-Normal\t$nbam\n$sampid-Tumor\t$tbam\n";  
		close OUT;

		open OUT,">$mydir/GATKCNV/config/somatic.list" or die "$!";;    
		print OUT "$sampid-Normal\t$sampid-Tumor\n";
		close OUT;
		my $cmd = "source /home/chengyuanfang/.bashrc && \\\n";
		$cmd.= "perl $configs{gatkcnvtool} -bam $mydir/GATKCNV/config/bamlst -sample_pair";
		$cmd.= " $mydir/GATKCNV/config/somatic.list -region $configs{gatkcnv_region}";
		$cmd.= " -outdir $mydir/GATKCNV/result $configs{gatkcnvparameter} && \\\n";
		my $cnvPNAME=$sampid."sCNV";
		$cmd.= "$configs{monitor} taskmonitor -q st.q -P $configs{priority} -p $cnvPNAME";
		$cmd.= " -i $mydir/GATKCNV/result/list/SomaticCNV_GATKACNV_dependence.txt";
		generateShell($shell,$cmd);
		$scripts{gatkcnv}{$sampid} = $shell;
	}
}

sub gatkcnv1{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
	my $mydir = "$outdir/$sampid/4.SCNV/gatkcnv";
	my $tbam = $hash{$sampid}{Tumor};
	my $nbam = $hash{$sampid}{Normal};
	#step1 CollectCounts of all samples 
	my $shell = "$outdir/$sampid/0.shell/a19.gatkcnv/s1.Normal.collectCounts.sh";
	`mkdir -p $outdir/$sampid/0.shell/a19.gatkcnv/tmp`;
	`mkdir -p $mydir/CollectCounts`;
	`mkdir -p $mydir/CNVcalling/$sampid`;
	my $cmd = "$configs{GATK4} --java-options \"-Xmx15G -Djava.io.tmpdir=$outdir/$sampid/0.shell/a19.gatkcnv/tmp\"";
	$cmd.= " CollectReadCounts -I $nbam -L $configs{gatkcnv_bed} -R $configs{reference}";
	$cmd.= " --format HDF5 --interval-merging-rule OVERLAPPING_ONLY";
	$cmd.= " --output $mydir/CollectCounts/$sampid-Normal.hdf5 && \\\n";
	$cmd.= "$configs{GATK4}  --java-options \"-Xmx15G -Djava.io.tmpdir=$outdir/$sampid/0.shell/a19.gatkcnv/tmp\"";
	$cmd.= " CollectAllelicCounts -I $nbam -L $configs{gatkcnv_bed}  -R $configs{reference}";
	$cmd.= " -O $mydir/CollectCounts/$sampid-Normal.allelicCounts.tsv";
	generateShell($shell,$cmd);
	$scripts{gatkcnv}{$sampid}{Normal} = $shell;
	$shell = "$outdir/$sampid/0.shell/a19.gatkcnv/s1.Tumor.collectCounts.sh";
	$cmd = "$configs{GATK4} --java-options \"-Xmx15G -Djava.io.tmpdir=$outdir/$sampid/0.shell/a19.gatkcnv/tmp\"";
	$cmd.= " CollectReadCounts -I $tbam -L $configs{gatkcnv_bed} -R $configs{reference}";
	$cmd.= " --format HDF5 --interval-merging-rule OVERLAPPING_ONLY --output $mydir/CollectCounts/$sampid-Tumor.hdf5 && \\\n";
	$cmd.= "$configs{GATK4} --java-options \"-Xmx15G -Djava.io.tmpdir=$outdir/$sampid/0.shell/a19.gatkcnv/tmp\"";
	$cmd.= " CollectAllelicCounts -I $tbam -L $configs{gatkcnv_bed}  -R $configs{reference}";
	$cmd.= " -O $mydir/CollectCounts/$sampid-Tumor.allelicCounts.tsv";
	generateShell($shell,$cmd);
	$scripts{gatkcnv}{$sampid}{Tumor} = $shell;
	#step2 CreateReadCountPanelOfNormals of normal samples
#   $shell = "$outdir/$sampid/0.shell/a19.gatkcnv/s2.CreateReadCountPanelOfNormals.sh";
	`ln -s $configs{pre_gatkcnvPoN} $mydir/cnv.pon.hdf5`;
#   generateShell($shell,$cmd);
#   $scripts{gatkcnv}{$sampid}{step2} = $shell;
	#step3 CNVcalling
	$shell = "$outdir/$sampid/0.shell/a19.gatkcnv/s3.run_cnv.sh";
	$cmd = "$configs{GATK4} --java-options \"-Xmx20g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
	$cmd.= " -XX:-UseGCOverheadLimit\" DenoiseReadCounts -I $mydir/CollectCounts/$sampid-Tumor.hdf5";
	$cmd.= " --count-panel-of-normals $mydir/cnv.pon.hdf5";
	$cmd.= " --standardized-copy-ratios $mydir/CNVcalling/$sampid-Tumor.standardizedCR.tsv";
	$cmd.= " --denoised-copy-ratios $mydir/CNVcalling/$sampid-Tumor.denoisedCR.tsv && \\\n";
	$cmd.= "$configs{GATK4} --java-options \"-Xmx100g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
	$cmd.= " -XX:-UseGCOverheadLimit\" ModelSegments --denoised-copy-ratios $mydir/CNVcalling/$sampid-Tumor.denoisedCR.tsv";
	$cmd.= " --allelic-counts $mydir/CollectCounts/$sampid-Tumor.allelicCounts.tsv";
	$cmd.= " --normal-allelic-counts $mydir/CollectCounts/$sampid-Normal.allelicCounts.tsv";
	$cmd.= " --output $mydir/CNVcalling/$sampid --output-prefix $sampid && \\\n";
	$cmd.= "$configs{GATK4} --java-options \"-Xmx20g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
	$cmd.= " -XX:-UseGCOverheadLimit\" CallCopyRatioSegments -I $mydir/CNVcalling/$sampid/$sampid.cr.seg";
	$cmd.= " -O $mydir/CNVcalling/$sampid/$sampid.called.seg && \\\n";
	$cmd.= "$configs{seg2vcf} $mydir/CNVcalling/$sampid/$sampid.called.seg $mydir/CNVcalling/$sampid/$sampid.called.vcf && \\\n";
	$cmd.= "perl $configs{gatkcnv_vep} --species homo_sapiens --merged --assembly GRCh37 --offline";
	$cmd.= " --pick_order canonical,tsl,biotype,rank,ccds,length --cache";
	$cmd.= " --dir /ldfssz1/ST_PRECISION/PMO/P18Z10200N0204/pd1/pipeline/chench/bin/software/VEP/v94/cache/";
	$cmd.= " --fasta /hwfssz4/BC_PUB/pipeline/DNA/DNA_Human_WES/DNA_Human_WES_2016b/Database/hg19/fa/hg19.fasta";
	$cmd.= " --format vcf --input_file $mydir/CNVcalling/$sampid/$sampid.called.vcf";
	$cmd.= " --output_file $mydir/CNVcalling/$sampid/$sampid.called.tab --tab --fork 4 --force_overwrite";
	generateShell($shell,$cmd);
	$scripts{gatkcnv}{$sampid}{step3} = $shell;
	#step4 ploting
	$shell = "$outdir/$sampid/0.shell/a19.gatkcnv/s4.plot.sh";
	$cmd = "$configs{GATK4} --java-options \"-Xmx2g\" PlotDenoisedCopyRatios";
	$cmd.= " --standardized-copy-ratios $mydir/CNVcalling/$sampid-Tumor.standardizedCR.tsv";
	$cmd.= " --denoised-copy-ratios $mydir/CNVcalling/$sampid-Tumor.denoisedCR.tsv";
	$cmd.= " --sequence-dictionary $configs{reference}.dict --minimum-contig-length 48129895";
	$cmd.= " --output $mydir/CNVcalling/$sampid --output-prefix $sampid && \\\n";
	$cmd.= "$configs{GATK4} --java-options \"-Xmx2g\" PlotModeledSegments";
	$cmd.= " --denoised-copy-ratios $mydir/CNVcalling/$sampid-Tumor.denoisedCR.tsv";
	$cmd.= " --allelic-counts $mydir/CNVcalling/$sampid/$sampid.hets.tsv";
	$cmd.= " --segments $mydir/CNVcalling/$sampid/$sampid.modelFinal.seg";
	$cmd.= " --sequence-dictionary $configs{reference}.dict --minimum-contig-length 48129895";
	$cmd.= " --output $mydir/CNVcalling/$sampid --output-prefix $sampid";
	generateShell($shell,$cmd);
	$scripts{gatkcnv}{$sampid}{step4} = $shell;

	}
}

sub integrate{
	my ($hash1,$hash2) = @_;
	my %hash1 = %{$hash1};
	my %hash2 = %{$hash2};
	for my $sampid(sort keys %hash1){
		my $mydir = "$outdir/$sampid/6.fusion/integrate";
		my $trnabam = $hash1{$sampid};
		my $tbam = $hash2{$sampid}{Tumor};
		my $nbam = $hash2{$sampid}{Normal};
		my $shell = "$outdir/$sampid/0.shell/a10.integrate/$sampid.integrate.sh";
		my $cmd = "$configs{INTEGRATE} fusion $configs{integratepar}";
		$cmd.= " -reads $mydir/reads.txt";
		$cmd.= " -sum $mydir/summary.tsv";
		$cmd.= " -ex $mydir/exons.tsv";
		$cmd.= " -bk $mydir/breakpoints.tsv";
		$cmd.= " -vcf $mydir/bk_sv.vcf";
		$cmd.= " -bedpe $mydir/fusions.bedpe";
		$cmd.= " -sample $sampid";
		$cmd.= " $configs{reference} $configs{INTEGRATEannot} $configs{INTEGRATEbwt} $trnabam $trnabam $tbam $nbam";
		generateShell($shell,$cmd);
		$scripts{integrate}{$sampid} = $shell;
		$fusionresult{$sampid}{integrate} = "$mydir/fusions.bedpe";
	}
}

sub starfusion{
	for my $sampid(sort keys %cleanfastq){
		my $mydir = "$outdir/$sampid/6.fusion/starfusion";
		open LIST, ">$mydir/fq.list" or die "$!";
		my @fq1 = @{$cleanfastq{$sampid}{fq1}};
		my @fq2 = @{$cleanfastq{$sampid}{fq2}};
		for my $i(0..$#fq1){
			print LIST join("\t",$sampid,$fq1[$i],$fq2[$i]),"\n";
		}
		close LIST;
		
		my $shell = "$outdir/$sampid/0.shell/a11.starfusion/$sampid.starfs.sh";
		my $cmd = "rm -rf $mydir/tmp\n";
		$cmd.= "$configs{STARfusion} --samples_file $mydir/fq.list --genome_lib_dir $configs{STARfslib}";
		$cmd.= " --CPU 6 --output_dir $mydir --bbmerge --STAR_SortedByCoordinate --STAR_PATH $configs{STAR}";
		$cmd.= " --outTmpDir $mydir/tmp --examine_coding_effect --extract_fusion_reads $configs{starfusionpar}";
		generateShell($shell,$cmd);
		$scripts{starfusion}{$sampid} = $shell;
		`rm -rf $mydir/tmp` if (-e "$mydir/tmp");
		$fusionresult{$sampid}{starfusion} = "$mydir/star-fusion.fusion_predictions.tsv";
	}
}

sub rsemstar{
	for my $sampid(sort keys %cleanfastq){
		my $mydir = "$outdir/$sampid/7.expression/RSEMstar";
		my $fq1 = join(",",@{$cleanfastq{$sampid}{fq1}});
		my $fq2 = join(",",@{$cleanfastq{$sampid}{fq2}});
		my $shell = "$outdir/$sampid/0.shell/a13.rsem/$sampid.rsem_star.sh";
		my $cmd = "$configs{RSEM} -p 2 --paired-end";
		$cmd.= " --phred33-quals" if ($phred{$sampid}{RNA} == 33);
		$cmd.= " --phred64-quals" if ($phred{$sampid}{RNA} == 64);
		$cmd.= " --star --star-path $configs{STARpath} $configs{rsempar}";
		$cmd.= " $fq1 $fq2 $configs{starref} $mydir/$sampid";
		generateShell($shell,$cmd);
		$scripts{rsem}{$sampid} = $shell;
		$exprresult{$sampid}{rsem} = "$mydir/$sampid.isoforms.results";
	}
}

sub readHLAlist{
	open LIST, "<$configs{hlalist}" or die "$!";
	while(my $line = <LIST>){
		chomp $line;
		my ($sampid,$hlafile) = (split /\t/,$line)[0,1];
		$hlatype{$sampid} = $hlafile;
	}
	close LIST;
}

sub hlasomatic{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/hlasomatic";
		for my $dtype(keys %{$hash{$sampid}}){
			if($dtype =~ /RNA/i){ next; }
			system("mkdir -p $mydir/$dtype/Polysolver") unless(-e "$mydir/$dtype/Polysolver");
			my $bam = $hash{$sampid}{$dtype};
			my $shell = "$outdir/$sampid/0.shell/a18.hlasomatic/$sampid.$dtype.Polysolver.sh";
			my $cmd = "source $configs{polysolver_config} && \\\n";
			$cmd.= "$configs{polysolver} $bam $configs{hlapar} $mydir/$dtype/Polysolver && \\\n";
			$cmd.= "perl $configs{combinescript} $mydir/$dtype/Polysolver $mydir/$sampid.$dtype.Polysolver.csv && \\\n";
			$cmd.= "rm -rf $mydir/$dtype/Polysolver/*.bam $mydir/$dtype/Polysolver/*.lik* $mydir/$dtype/Polysolver/*.sam";
			generateShell($shell, $cmd);
			$scripts{polysolver}{$sampid}{$dtype} = $shell;
			$hlatype{$sampid}{$dtype} = "$mydir/$sampid.$dtype.Polysolver.csv";
		}
	}
}

sub HLAminer{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			if ($dtype =~ /RNA/i){ next; }
			my $mydir = "$outdir/$sampid/5.somatic/hlasomatic/$dtype/HLAminer";
			system("mkdir -p $mydir") unless (-e $mydir);
			open FOF, "> $mydir/$sampid.$dtype.fof";
			# print FOF "$outdir/$sampid/2.sentieon/$dtype/$sampid.$dtype.HLA.fq.gz\n";
			print FOF "$outdir/$sampid/2.megaBOLT/$sampid-$dtype/HLA.fq.gz\n";
			print FOF "$unmapfqs{$sampid}{$dtype}\n";
			close FOF;
			my $shell = "$outdir/$sampid/0.shell/a18.hlasomatic/$sampid.$dtype.HLAminer.sh";
			my $cmd = "###TASR\necho \"Running TASR...\"\n";
			$cmd.= "$configs{HLAminer}/TASR -f $mydir/$sampid.$dtype.fof -s $configs{TASRtargets}";
			$cmd.= " -b $mydir/$sampid.$dtype.TASRhla $configs{TASRpar} && \\\n";
			$cmd.= "###Restrict 200nt+ contigs\n";
			$cmd.= "cat $mydir/$sampid.$dtype.TASRhla.contigs |";
			$cmd.= " perl -ne 'if(/size(\\d+)/){if(\$1>=200){\$flag=1;print;}else{\$flag=0;}}else{print if(\$flag);}' >";
			$cmd.= " $mydir/$sampid.$dtype.TASRhla200.contigs && \\\n";
			$cmd.= "###Create a [NCBI] blastable database\n";
			$cmd.= "echo \"Formatting blastable database...\"\n";
			$cmd.= "$configs{HLAminer}/formatdb -p F -i $mydir/$sampid.$dtype.TASRhla200.contigs -l $mydir/$sampid.$dtype.formatdb.log && \\\n";
			$cmd.= "###Align contigs against database\n";
			$cmd.= "echo \"Aligning TASR contigs to HLA references...\"\n";
			$cmd.= "$configs{HLAminer}/parseXMLblast.pl -c $configs{HLAminer}/ncbiBlastConfig.txt -d $configs{TASRtargets}";
			$cmd.= " -i $mydir/$sampid.$dtype.TASRhla200.contigs -o 0 -a 1 > $mydir/$sampid.$dtype.tig_vs_hla-ncbi.coord && \\\n";
			$cmd.= "###Align HLA references to contigs\n";
			$cmd.= "echo \"Aligning HLA references to TASR contigs (go have a coffee, it may take a while)...\"\n";
			$cmd.= "$configs{HLAminer}/parseXMLblast.pl -c $configs{HLAminer}/ncbiBlastConfig.txt -d $mydir/$sampid.$dtype.TASRhla200.contigs";
			$cmd.= " -i $configs{TASRtargets} -o 0 -a 1 > $mydir/$sampid.$dtype.hla_vs_tig-ncbi.coord && \\\n";
			$cmd.= "###Predict HLA alleles\n";
			$cmd.= "echo \"Predicting HLA alleles...\"\n";
			$cmd.= "$configs{HLAminer}/HLAminer.pl -b $mydir/$sampid.$dtype.tig_vs_hla-ncbi.coord -r $mydir/$sampid.$dtype.hla_vs_tig-ncbi.coord";
			$cmd.= " -c $mydir/$sampid.$dtype.TASRhla200.contigs -h $configs{TASRtargets} -m $mydir -p $configs{hlanomp} && \\\n";
			$cmd.= "$configs{HLAminer_Parse} $mydir/HLAminer_HPTASR.csv > $outdir/$sampid/5.somatic/hlasomatic/$sampid.$dtype.HLAminer_HPTASR.csv";
			generateShell($shell, $cmd);
			$scripts{HLAminer}{$sampid}{$dtype} = $shell;
		}
	}
}

sub NeoantigenPredict_rna{
	for my $sampid(sort keys %snvvcf){
		my $snvvcf = (exists $snvvcf{$sampid}{merge})? $snvvcf{$sampid}{merge} : $snvvcf{$sampid}{mutect};
		my $indelvcf = (exists $indelvcf{$sampid}{merge})? $indelvcf{$sampid}{merge} : $indelvcf{$sampid}{strelka2};
		my $fusionbedpe = $fusionresult{$sampid}{integrate};
		my $exprresult = $exprresult{$sampid}{rsem};
		my $hla = $hlatype{$sampid};
		my $mydir = "$outdir/$sampid/8.neoantigen";
		## SNV
		my $shell = "$outdir/$sampid/0.shell/a15.neoantigen/run.snv.sh";
		my $cmd = "less $snvvcf | cut -f1,2,3,4,5,6,7 > $outdir/$sampid/8.neoantigen/snv.mut.vcf && \\\n";
		$cmd.= "$configs{python3} $configs{Neoantigenpy} -a $hla -m $configs{methodlist}";
		$cmd.= " -t snv -e $configs{assembly} -o $mydir -min $configs{minLength} -max $configs{maxLength}";
		$cmd.= " -i $outdir/$sampid/8.neoantigen/snv.mut.vcf $configs{nt} -s snv -exp $exprresult";
		generateShell($shell,$cmd);
		$scripts{neoantigen}{$sampid}{snv} = $shell;
		## InDel
		$shell = "$outdir/$sampid/0.shell/a15.neoantigen/run.indel.sh";
		$cmd = "less $indelvcf | cut -f1,2,3,4,5,6,7 > $outdir/$sampid/8.neoantigen/indel.mut.vcf && \\\n";
		$cmd.= "$configs{python3} $configs{Neoantigenpy} -a $hla -m $configs{methodlist} -t indel";
		$cmd.= " -e $configs{assembly} -o $mydir -min $configs{minLength} -max $configs{maxLength}";
		$cmd.= " -i $outdir/$sampid/8.neoantigen/indel.mut.vcf $configs{nt} -s indel -exp $exprresult";
		generateShell($shell,$cmd);
		$scripts{neoantigen}{$sampid}{indel} = $shell;
		## Fusion
		if ($configs{integrate} =~ /true/i){
			$shell = "$outdir/$sampid/0.shell/a15.neoantigen/run.fusion.sh";
			$cmd = "$configs{python3} $configs{Neoantigenpy} -a $hla -m $configs{methodlist}";
			$cmd.= " -t fusion -e $configs{assembly} -o $mydir -min $configs{minLength} -max $configs{maxLength}";
			$cmd.= " -i $fusionbedpe $configs{nt} -s fusion -exp $exprresult";
			generateShell($shell,$cmd);
			$scripts{neoantigen}{$sampid}{fusion} = $shell;
		}
	}
}

sub NeoantigenPredict{
	for my $sampid(sort keys %snvvcf){
		my $snvvcf = (exists $snvvcf{$sampid}{merge})? $snvvcf{$sampid}{merge} : $snvvcf{$sampid}{mutect};
		my $indelvcf = (exists $indelvcf{$sampid}{merge})? $indelvcf{$sampid}{merge} : $indelvcf{$sampid}{strelka2};
		my $fusionbedpe = $fusionresult{$sampid}{integrate};
		my $hla = $hlatype{$sampid};
		my $mydir = "$outdir/$sampid/8.neoantigen";
		## SNV
		my $shell = "$outdir/$sampid/0.shell/a15.neoantigen/run.snv.sh";
		my $cmd = "less $snvvcf | cut -f1,2,3,4,5,6,7 > $outdir/$sampid/8.neoantigen/snv.mut.vcf && \\\n";
		$cmd.= "$configs{python3} $configs{Neoantigenpy} -a $hla -m $configs{methodlist}";
		$cmd.= " -t snv -e $configs{assembly} -o $mydir -min $configs{minLength} -max $configs{maxLength}";
		$cmd.= " -i $outdir/$sampid/8.neoantigen/snv.mut.vcf $configs{nt} -s snv";
		generateShell($shell,$cmd);
		$scripts{neoantigen}{$sampid}{snv} = $shell;
		## InDel
		$shell = "$outdir/$sampid/0.shell/a15.neoantigen/run.indel.sh";
		$cmd = "less $indelvcf | cut -f1,2,3,4,5,6,7 > $outdir/$sampid/8.neoantigen/indel.mut.vcf && \\\n";
		$cmd.= "$configs{python3} $configs{Neoantigenpy} -a $hla -m $configs{methodlist} -t indel";
		$cmd.= " -e $configs{assembly} -o $mydir -min $configs{minLength}";
		$cmd.= " -max $configs{maxLength} -i $outdir/$sampid/8.neoantigen/indel.mut.vcf $configs{nt} -s indel";
		generateShell($shell,$cmd);
		$scripts{neoantigen}{$sampid}{indel} = $shell;
		## Fusion
		if($configs{integrate} =~ /true/){
			$shell = "$outdir/$sampid/0.shell/a15.neoantigen/run.fusion.sh";
			$cmd = "$configs{python3} $configs{Neoantigenpy} -a $hla -m $configs{methodlist}";
			$cmd.= " -t fusion -e $configs{assembly} -o $mydir -min $configs{minLength}";
			$cmd.= " -max $configs{maxLength} -i $fusionbedpe $configs{nt} -s fusion";
			generateShell($shell,$cmd);
			$scripts{neoantigen}{$sampid}{fusion} = $shell;
		}
	}
}

sub chemicalDrug{
	my (%hash) = @_;
	for my $sampid(sort keys %hash){
		my $tbam = $hash{$sampid}{Tumor};
		my $nbam = $hash{$sampid}{Normal};
		for my $dtype(keys %{$hash{$sampid}}){
			if($dtype =~ /Normal/i){
				my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/a02.annosnp_$dtype.sh";
				my $cmd = "mkdir -p $outdir/$sampid/2.megaBOLT/$dtype && \\\n";
				$cmd.= "java -Xmx8G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$outdir/gatksnp/tmp";
				$cmd.= " -jar $configs{GATK3} -T HaplotypeCaller";
				$cmd.= " -R $configs{reference} -I $nbam $configs{happar}";
				$cmd.= " -o $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.sortdup.bqsr.bam.HaplotypeCaller.vcf && \\\n";
				$cmd.= "$configs{java} -Djava.io.tmpdir=$outdir/$sampid/2.megaBOLT/$sampid-$dtype/tmp";
				$cmd.= " -jar $configs{GATK3} -T VariantFiltration -R ";
				$cmd.= "$configs{reference} $configs{SNPfilterParameter} ";
				$cmd.= "$outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.sortdup.bqsr.bam.HaplotypeCaller.vcf";
				$cmd.= " --logging_level ERROR -o $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.sortdup.bqsr.vcf && \\\n";
				$cmd.= "awk '\$1~/^#/ || \$7==\"PASS\"' $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.sortdup.bqsr.vcf";
				$cmd.= " > $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.final.variants.vcf && \\\n";
				$cmd.= "perl $configs{ANNOVAR} $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.final.variants.vcf";
				$cmd.= " $configs{ANNOVAR_refdir} -out $outdir/$sampid/2.megaBOLT/$sampid-$dtype/final.variants $configs{annovar_par} && \\\n";
				$cmd.= "less $outdir/$sampid/2.megaBOLT/$sampid-$dtype/final.variants.hg19_multianno.txt | cut -f1,2,4,5,6 |";
				$cmd.= " grep 'rs' > $outdir/$sampid/2.megaBOLT/$sampid-$dtype/final.variants.txt";
				generateShell($shell,$cmd);
				$scripts{gatksnp_ann}{$sampid}{$dtype}=$shell;
				##snp chemicalDrug
				$shell = "$outdir/$sampid/0.shell/a16.medicine/a16.chemicalDrug_$dtype.sh";
				$cmd = "perl $configs{chemicaldrug_tool} $outdir/$sampid/2.megaBOLT/$sampid-$dtype/final.variants.txt";
				$cmd.= " > $outdir/$sampid/9.medicine/chemicalDrug/chemicalDrug.out";
				generateShell($shell,$cmd);
				$scripts{chemicalDurg_ann}{$sampid}{$dtype}=$shell;
			}
		}
	}
}
sub resistantDrug{
	for my $sampid(sort keys %snvvcf){
		my $shell="$outdir/$sampid/0.shell/a16.medicine/resistantDrug.sh";
		my $cmd = "perl $configs{resistant_tool} $configs{purity} $snvvcf{$sampid}{target_maf}";
		$cmd.= " $outdir/$sampid/9.medicine/resistantDrug && \\\n";
		$cmd.= "cat $outdir/$sampid/5.somatic/result/softmerge.snv.vcf";
		$cmd.= " $outdir/$sampid/5.somatic/result/softmerge.indel.vcf";
		$cmd.= " >$outdir/$sampid/5.somatic/result/softmerge.snv_indel.vcf && \\\n";
		$cmd.= "perl $configs{resistant_freq} $outdir/$sampid/9.medicine/resistantDrug/cancer_gene.raw.txt";
		$cmd.= " $outdir/$sampid/5.somatic/result/softmerge.snv_indel.vcf";
		generateShell($shell,$cmd);
		$scripts{resist}{$sampid} = $shell;
	}
}

sub immupathway{
	for my $sampid(sort keys %snvvcf){
		my $shell="$outdir/$sampid/0.shell/a17.immupathway/run.sh";
		my $mydir="$outdir/$sampid/10.immupathway";
		my $cmd ="perl $configs{immupathway_tool1} $snvvcf{$sampid}{target_maf}";
		$cmd.= " $mydir/immunepathway $sampid $configs{purity} \n";
		$cmd.= "perl $configs{immupathway_tool2} $mydir/immunepathway/$sampid.SNV_SIndel.refGene.list";
		$cmd.= " > $mydir/immunepathway/$sampid.SNV_SIndel.refGene.glist\n";
		$cmd.= "$configs{Rscript} $configs{immupathway_tool3} $mydir/immunepathway/$sampid.SNV_SIndel.refGene.glist";
		$cmd.= " $mydir/immunepathway/$sampid.SNV_SIndel.refGene.glist.Rscript\n";
		$cmd.= "ls $mydir/immunepathway/*_enrichKEGG.txt |";
		$cmd.= " perl -ne 'chomp;\$name=(split /\\//)[-1];\$name=~s/(.*).SNV_SIndel.refGene.glist_enrichKEGG.txt/\$1/;print";
		$cmd.= " \"\$name\\t\$_\\n\";' > $mydir/immunepathway/enrichKEGG.list\n";
		$cmd.= "perl $configs{immupathway_tool4} $mydir/immunepathway/enrichKEGG.list $mydir/immunepathway/enrichKEGG.txt\n";
		$cmd.= "perl $configs{immupathway_tool5} $mydir/immunepathway $configs{immupathway_data}";
		$cmd.= " > $mydir/immunepathway/$sampid.enrichKEGG.final.out";
		generateShell($shell,$cmd);
		$scripts{immunepathway}{$sampid} = $shell;
	}
}

sub getreport{
	for my $sampid(sort keys %snvvcf){
		my $shell="$outdir/$sampid/0.shell/step1.getreportinfo.sh";
		my $cmd ="perl $configs{getreport_too11} $outdir/$sampid >$outdir/$sampid/$sampid.reportinfo.xls";
		generateShell($shell,$cmd);
		$scripts{getreport}{$sampid} = $shell;
	}
}
	
sub edgeList{
	my ($day) = @_;
	my $mydir = "$outdir/shell_run";
	open EDGE, ">$outdir/shell_run/edge.$configs{projectname}.$day.list" or die $!;
	%fqinputs = %bqsrbam if ($configs{alignment} =~ /false/i);
	for my $sampid(sort keys %fqinputs){
		for my $dtype(sort keys %{$fqinputs{$sampid}}){
			if (exists $scripts{clean}){        
				for my $lane(sort keys %{$scripts{clean}{$sampid}{$dtype}}){
					foreach my $num(sort keys %{$scripts{clean}{$sampid}{$dtype}{$lane}}){
						if($configs{split_num} > 1 && $splitfqs{$sampid}{$dtype}{$lane}{fqstat_exist} =~ /true/i){
							## split fq
							my $t = $num % $configs{split_num};
							my $zui = ($num - $t) / $configs{split_num}; 
							print EDGE "$scripts{split}{$sampid}{$dtype}{$lane}{$zui.1}:1G:1CPU $scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU\n";
							print EDGE "$scripts{split}{$sampid}{$dtype}{$lane}{$zui.2}:1G:1CPU $scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU\n";   
						}
						## clean megaBOLT
						print EDGE "$scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU $scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q}\n"; 
					}
			   }
			}

			## alignment and qc
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{megaboltqc}{$sampid}:10G:1CPU\n" if($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /true/i && $configs{tumor_only} =~ /false/i && -e $configs{qcvcf});
			print EDGE "$scripts{megaboltqc}{$sampid}:10G:1CPU\n" if ($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /false/i && $dtype =~ /Normal/ && $configs{tumor_only} =~ /false/i && -e $configs{qcvcf});
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{coverage}{$sampid}{$dtype}:10G:1CPU\n" if($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{coverage}{$sampid}{$dtype}:10G:1CPU\n" if ($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /false/i);
			## microbiome analysis
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{metaphlan}{$sampid}{$dtype}:$configs{metaphlan_q}\n" if($configs{metaphlan} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{metaphlan}{$sampid}{$dtype}:$configs{metaphlan_q}\n" if($configs{metaphlan} =~ /true/i && $configs{alignment} =~ /false/i);
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{kraken2}{$sampid}{$dtype}:$configs{kraken2_q}\n" if($configs{kraken} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{kraken2}{$sampid}{$dtype}:$configs{kraken2_q}\n" if($configs{kraken} =~ /true/i && $configs{alignment} =~ /false/i);                
			## msi
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{mantis}{$sampid}:$configs{mantis_q}\n" if($configs{msi} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{mantis}{$sampid}:$configs{mantis_q}\n" if($configs{msi} =~ /true/i && $configs{alignment} =~ /false/i && $dtype =~ /Normal/);
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{msisensor}{$sampid}:10G:1CPU\n" if($configs{msi} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{msisensor}{$sampid}:10G:1CPU\n" if($configs{msi} =~ /true/i && $configs{alignment} =~ /false/i && $dtype =~ /Normal/);
			## gatk_snp_ann
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{gatksnp_ann}{$sampid}{$dtype}:1G:1CPU\n" if($configs{alignment} =~ /true/i);
			print EDGE "$scripts{gatksnp_ann}{$sampid}{$dtype}:1G:1CPU $scripts{chemicalDurg_ann}{$sampid}{$dtype}:1G:1CPU\n" if($configs{chemicalDrug} =~ /true/i);
			## SV
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{Manta}{$sampid}:45G:8CPU\n" if($configs{alignment} =~ /true/i && $configs{manta} =~ /true/i);
			print EDGE "$scripts{Manta}{$sampid}:45G:8CPU\n" if($configs{manta} =~ /true/i && $configs{alignment} =~ /false/i);

			for my $mychr(@chrs){
				## Mutect
				if($configs{mutect} =~ /true/i){
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{mutect}{step1}{$sampid}{$mychr}:10G:1CPU\n" if($configs{alignment} =~ /true/i);
					print EDGE "$scripts{mutect}{step1}{$sampid}{$mychr}:10G:1CPU $scripts{mutect}{step2}{$sampid}:1G:1CPU\n" if($dtype =~ /Normal/);
				}
				## varscan2
				if($configs{varscan2} =~ /true/i){
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{varscan2}{step1}{$sampid}{$mychr}:6G:1CPU\n" if($configs{alignment} =~ /true/i);
					print EDGE "$scripts{varscan2}{step1}{$sampid}{$mychr}:6G:1CPU $scripts{varscan2}{step2}{$sampid}:1G:1CPU\n" if($dtype =~ /Normal/);
				}
				## MuSE
				if($configs{MuSE} =~ /true/i){
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{MuSE}{step1}{$sampid}{$mychr}:8G:1CPU\n" if($configs{alignment} =~ /true/i);
					print EDGE "$scripts{MuSE}{step1}{$sampid}{$mychr}:8G:1CPU $scripts{MuSE}{step2}{$sampid}:20G:1CPU\n" if($dtype =~ /Normal/);
				}
				## svaba
				if($configs{svaba} =~ /true/i){
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{svaba}{$sampid}{$mychr}:10G:1CPU\n" if($configs{alignment} =~ /true/i);
					print EDGE "$scripts{svaba}{$sampid}{$mychr}:10G:1CPU $scripts{svaba}{$sampid}{merge}:1G:1CPU\n" if ($dtype =~ /Normal/);
				}
				## strelka
				if($configs{strelka} =~ /true/i && $dtype =~ /Normal/){
					print EDGE "$scripts{strelka}{step1}{$sampid}:8G:1CPU $scripts{strelka}{step2}{$sampid}{$mychr}:8G:1CPU\n";
					print EDGE "$scripts{strelka}{step2}{$sampid}{$mychr}:8G:1CPU $scripts{strelka}{step3}{$sampid}:1G:1CPU\n";
				}
				## mutect2               
				if($configs{mutect2} =~ /true/i){
					print EDGE "$scripts{mutect2}{$sampid}{PoN}:15G:2CPU $scripts{mutect2}{step1}{$sampid}{$mychr}:8G:1CPU\n" unless($configs{panel_of_normal});
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{mutect2}{step1}{$sampid}{$mychr}:8G:1CPU\n" if($configs{alignment} =~ /true/i && $configs{panel_of_normal});
					print EDGE "$scripts{mutect2}{step1}{$sampid}{$mychr}:8G:1CPU $scripts{mutect2}{step2}{$sampid}:1G:1CPU\n" if($dtype =~ /Normal/);
				}
			}
				
			## Strelka2 (BQSR bam)
			if($configs{strelka2} =~ /true/i){
				print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{strelka2}{$sampid}:6G:4CPU\n" if ($configs{alignment} =~ /true/i);
				print EDGE "$scripts{strelka2}{$sampid}:6G:4CPU\n" if ($configs{alignment} =~ /false/i && $dtype =~ /Normal/);
			}
			## Mutect2 (BQSR bam)
			if($configs{mutect2} =~ /true/i){
				print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{mutect2}{$sampid}{contamination}:6G:1CPU\n" if($configs{alignment} =~ /true/i);
				print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{mutect2}{$sampid}{artifact}:4G:1CPU\n" if ($configs{alignment} =~ /true/i && $dtype =~ /Tumor/i);
				print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{mutect2}{$sampid}{PoN}:15G:2CPU\n" if($configs{alignment} =~ /true/i && not $configs{panel_of_normal});
				# print EDGE "$scripts{mutect2}{$sampid}{contamination}:6G:1CPU\n" if($configs{alignment} =~ /false/i && $dtype =~ /Normal/);
				# print EDGE "$scripts{mutect2}{$sampid}{artifact}:4G:1CPU\n"if ($configs{alignment} =~ /false/i && $dtype =~ /Tumor/i);
				# print EDGE "$scripts{mutect2}{$sampid}{PoN}:15G:2CPU\n" if ($configs{alignment} =~ /false/i && $dtype =~ /Normal/);
				print EDGE "$scripts{mutect2}{$sampid}{contamination}:6G:1CPU $scripts{mutect2}{$sampid}{filter}:25G:1CPU\n" if($dtype =~ /Normal/);
				print EDGE "$scripts{mutect2}{$sampid}{artifact}:4G:1CPU $scripts{mutect2}{$sampid}{filter}:25G:1CPU\n" if($dtype =~ /Normal/);
				print EDGE "$scripts{mutect2}{step2}{$sampid}:1G:1CPU $scripts{mutect2}{$sampid}{filter}:25G:1CPU\n" if($dtype=~/Normal/);
			}
			## Strelka (BQSR bam)
			if($configs{strelka} =~ /true/i){
				print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{strelka}{step1}{$sampid}:8G:1CPU\n" if($configs{alignment} =~ /true/i);
			}
			## tumor_only result postprocessing
			if($configs{tumor_only} =~ /true/i && $dtype =~ /Tumor/i){
				print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{tumor_only}{$sampid}{contamination}:6G:1CPU\n" if($configs{alignment} =~ /true/i);
				# print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{tumor_only}{$sampid}{artifact}:10G:1CPU\n" if ($configs{alignment} =~ /true/i);
				print EDGE "$scripts{tumor_only}{$sampid}{contamination}:6G:1CPU $scripts{tumor_only}{$sampid}{filter}:16G:1CPU\n";
				# print EDGE "$scripts{tumor_only}{$sampid}{artifact}:10G:1CPU $scripts{tumor_only}{$sampid}{filter}:16G:1CPU\n";
			}
			## scnv
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{FACETS_pileup}{$sampid}:3G:1CPU\n" if ($configs{FACETS} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{FACETS_pileup}{$sampid}:3G:1CPU $scripts{FACETS_call}{$sampid}:3G:1CPU\n" if($configs{FACETS}=~/true/i && $dtype=~/Normal/);
			print EDGE  "$scripts{FACETS_call}{$sampid}:3G:1CPU $scripts{FACETS_anno}{$sampid}:3G:1CPU\n" if($configs{FACETS}=~/true/i && $dtype=~/Normal/); 
			# print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{gatkcnv}{$sampid}:1G:1CPU\n" if($configs{gatkcnv} =~ /true/i && $configs{alignment} =~ /true/i);
			# print EDGE "$scripts{gatkcnv}{$sampid}:1G:1CPU\n" if($configs{gatkcnv} =~ /true/i && $configs{alignment} =~ /false/i && $dtype =~ /Normal/);
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{gatkcnv}{$sampid}{$dtype}:15G:1CPU\n" if($configs{gatkcnv} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{gatkcnv}{$sampid}{$dtype}:15G:1CPU $scripts{gatkcnv}{$sampid}{step3}:100G:1CPU\n" if($configs{gatkcnv} =~ /true/i);
			print EDGE "$scripts{gatkcnv}{$sampid}{step3}:100G:1CPU $scripts{gatkcnv}{$sampid}{step4}:5G:1CPU\n" if($configs{gatkcnv} =~ /true/i && $dtype =~ /Normal/);
			## INTEGRATE (re-aligned bam)
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{integrate}{$sampid}:25G:1CPU\n" if ($configs{integrate} =~ /true/i && $configs{alignment} =~ /true/i); 
				
			## targetdrug
			if ($configs{NeoantigenPredict} =~ /true/i){
				print EDGE "$scripts{neoantigen}{$sampid}{snv}:5G:6CPU $scripts{targetDrug}{$sampid}:1G:1CPU\n" if($configs{targetDrug} =~ /true/i && $dtype =~ /Normal/);
			}
			print EDGE "$scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU $scripts{targetDrug}{$sampid}:1G:1CPU\n" if($configs{targetDrug} =~ /true/i && $dtype =~ /Normal/);                
			##immunepathway
			print EDGE "$scripts{targetDrug}{$sampid}:1G:1CPU $scripts{immunepathway}{$sampid}:1G:1CPU\n" if($configs{immupathway} =~ /true/i && $dtype =~ /Normal/);
			##hla dection
			if($configs{alignment} =~ /true/i){

			}
			if($configs{alignment} =~ /false/i){
				unless(-e $configs{hlalist}){
					## polysolver
					print EDGE "$scripts{polysolver}{$sampid}{$dtype}:4G:1CPU\n" if($configs{hlasomatic} =~ /true/i);
					## HLAminer
					print EDGE "$scripts{HLAminer}{$sampid}{$dtype}:$configs{HLAminer_q}\n" if($configs{hlaminer} =~ /true/i);
				}
			}
			else{
				unless(-e $configs{hlalist}){
					## polysolver
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{polysolver}{$sampid}{$dtype}:4G:1CPU\n" if($configs{hlasomatic} =~ /true/i);
					## HLAminer
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{HLAminer}{$sampid}{$dtype}:$configs{HLAminer_q}\n" if($configs{hlaminer} =~ /true/i);
					
				}
			}
			## Neoantigen Prediction
			# print EDGE "$scripts{mutect}{step2}{$sampid}:1G:1CPU $scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU\n" if ($configs{mutect} =~ /true/i && $configs{mutect2} =~ /true/i && $dtype =~ /Normal/);
			# print EDGE "$scripts{mutect2}{$sampid}{filter}:25G:1CPU $scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU\n" if ($configs{mutect} =~ /true/i && $configs{mutect2} =~ /true/i && $dtype =~ /Normal/);
			if ($configs{somatic} =~ /true/i && $dtype =~ /Normal/i){
				for my $mychr(@chrs){
					print EDGE "$scripts{mutect}{step2}{$sampid}:1G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{mutect} =~ /true/i && $configs{mutect2} =~ /true/i);
					print EDGE "$scripts{mutect2}{$sampid}{filter}:35G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{mutect} =~ /true/i && $configs{mutect2} =~ /true/i);
					print EDGE "$scripts{strelka}{step3}{$sampid}:6G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{strelka} =~ /true/i);
					print EDGE "$scripts{MuSE}{step2}{$sampid}:20G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{MuSE} =~ /true/i);
					print EDGE "$scripts{svaba}{$sampid}{merge}:1G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{svaba} =~/true/i);
					print EDGE "$scripts{strelka2}{$sampid}:6G:1CPU $scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU\n" if ($configs{strelka2}=~/true/i);
					print EDGE "$scripts{mergeSnvInDel}{$sampid}{$mychr}:5G:1CPU $scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU\n";
				}
				print EDGE "$scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU $scripts{phase}{$sampid}:2G:2CPU\n";
			}
			if ($configs{NeoantigenPredict} =~ /true/i){
				print EDGE "$scripts{hlanormal}{$sampid}:4G:1CPU $scripts{neoantigen}{$sampid}{snv}:5G:6CPU\n" if(!(-e $configs{hlalist}) && $dtype =~/Normal/);
				print EDGE "$scripts{hlanormal}{$sampid}:4G:1CPU $scripts{neoantigen}{$sampid}{indel}:5G:6CPU\n" if(!(-e $configs{hlalist}) && $dtype =~/Normal/);
				print EDGE "$scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU $scripts{neoantigen}{$sampid}{snv}:5G:6CPU\n" if($dtype =~/Normal/);
				print EDGE "$scripts{mergeSnvInDel}{$sampid}{concat}:5G:2CPU $scripts{neoantigen}{$sampid}{indel}:5G:6CPU\n" if($dtype =~/Normal/);
				if ($configs{integrate} =~ /true/i){
					print EDGE "$scripts{hlanormal}{$sampid}:4G:1CPU $scripts{neoantigen}{$sampid}{fusion}:1G:6CPU\n" unless(-e $configs{hlalist});
					print EDGE "$scripts{integrate}{$sampid}:25G:1CPU $scripts{neoantigen}{$sampid}{fusion}:1G:6CPU\n";
				}
			}
			## resist drug
			print EDGE "$scripts{targetDrug}{$sampid}:1G:1CPU $scripts{resist}{$sampid}:1G:1CPU\n" if($configs{resistant} =~ /true/i && $dtype =~/Normal/);
		}
	}
	close EDGE;
}