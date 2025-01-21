#! /usr/bin/env perl -w
# @Author: Xiangyu Guan
# @Date:   2025-01-21
# @Last Modified by:   Xiangyu Guan
# @Last Modified time: 2025-01-21
use strict;
use warnings;
use POSIX;
use File::Basename;
use Data::Dumper;

@ARGV == 3 or die "Usage: perl $0 <Fastq.list> <Configs> <Outdir>\n";
my ($vcflist,$config,$outdir) = @ARGV;
system("mkdir -p $outdir") unless(-e $outdir);
#######################################################################
my (%configs, %scripts);
my (%samples, %snpvcf, %snvvcf, %indelvcf, %Kmer, %CdsMutPro, %BAM);
my (@Blocks, @ns_len, @pep_len);

########################### Main Parameters ###########################
readConfigs($config);
readVcfList($vcflist);

open BLOCK, "<$configs{chrom_block}";
while (my $line=<BLOCK>) {
    next if($line =~ /^#/);
    chomp $line;
    push @Blocks, $line;
}
close BLOCK;
push @Blocks, 'other';
for (my $len = $configs{ns_min_length}; $len <= $configs{ns_max_length}; $len += 3) {
# for my $len($configs{ns_min_length}..$configs{ns_max_length}){
    # next unless($len % 3 == 0);
    push @ns_len, $len;
    push @pep_len, int($len / 3);
}
$configs{pep_length_range} = join(" ", @pep_len) unless($configs{pep_length_range});

CdsDatabase(%samples);
WGSDatabase(%samples);

# Generate day info
chomp(my $day = `date +%Y%m%d`);
## run scripts
EdgeList($day);
## Generate run.sh
open RUNSH, ">$outdir/shell_run/run.$configs{projectname}.$day.sh" or die "$!";
print RUNSH "$configs{monitor} taskmonitor --q1 $configs{queue} --q2 $configs{queue2} --P1 $configs{priority} --P2 $configs{priority2}";
print RUNSH " -p $configs{projectname} -i $outdir/shell_run/edge.$configs{projectname}.$day.list -f 1\n";
close RUNSH;
########################### Functions ###########################
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
    print OUT "qstat -j $shell_name > $output_shell.log\n";
    close OUT;
}

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

sub readVcfList{
    my ($inputs) = @_;
    open INPUT, "<$inputs" or die "$!";
    my $num = 0;
    while(my $line = <INPUT>){
        chomp $line;
        # my ($sampid,$snp,$snv,$indel, $kmer_path) = split /\s+/,$line;
        my ($sampid,$snp,$snv,$indel,$bam) = split /\s+/,$line;
        $samples{$sampid} = 1;
        $snpvcf{$sampid} = $snp unless($snp eq "-");
        $snvvcf{$sampid} = $snv unless($snv eq "-");
        $indelvcf{$sampid} = $indel unless($indel eq "-");
        # $Kmer{$sampid} = $kmer_path unless($kmer_path eq "-");
        $BAM{$sampid} = $bam unless($bam eq "-");
    }
    close INPUT;
}

sub CdsDatabase{
    my (%hash) = @_;
    for my $sampid(keys %hash){
        my $mydir = "$outdir/$sampid/1.neoantigen/CDS";
        system("mkdir -p $mydir") unless(-e $mydir);
        system("mkdir -p $outdir/$sampid/0.shell/a01.neoantigen") unless(-e "$outdir/$sampid/0.shell/a01.neoantigen");
        my $shell = "$outdir/$sampid/0.shell/a01.neoantigen/a01.$sampid.MutCDS.sh";
        my $cmd = "$configs{CDS_DB_pl}";
        $cmd.= " -s $snvvcf{$sampid}" if($snvvcf{$sampid});
        $cmd.= " -i $indelvcf{$sampid}" if($indelvcf{$sampid});
        $cmd.= " -g $snpvcf{$sampid}" if($snpvcf{$sampid});
        $cmd.= " -c $configs{cds_from_genomic}" if($configs{cds_from_genomic});
        $cmd.= " -a $configs{chromAlias}" if($configs{chromAlias});
        $cmd.= " -p $mydir/$sampid && \\\n";
        $cmd.= "$configs{CDS_Uniprot_pl}";
        $cmd.= " -r $configs{Protein_ref}" if($configs{Protein_ref});
        $cmd.= " -i $mydir/$sampid.cds.faa.gz -p $mydir/$sampid && \\\n";
        $cmd.= "$configs{CDS_Var2Tsv_pl} $configs{cds_from_genomic} $configs{chromAlias}";
        $cmd.= " $mydir/$sampid.cds.faa.gz $mydir/$sampid.MutCDS.tsv.gz";
        generateShell($shell, $cmd);
        $scripts{$sampid}{cds_db} = $shell;
        $CdsMutPro{$sampid}{cds_faa} = "$mydir/$sampid.cds.faa.gz";
        $CdsMutPro{$sampid}{cds_plus_uniprot_db} = "$mydir/$sampid.protein.DB.fasta.gz";
    }
}

sub WGSDatabase{
    my (%hash) = @_;
    for my $sampid(keys %hash){
        my $mydir = "$outdir/$sampid/1.neoantigen/WGS";
        system("mkdir -p $mydir") unless(-e $mydir);
        my $shell = "$outdir/$sampid/0.shell/a01.neoantigen/a01.$sampid.MutWGS.sh";
        my $cmd = "$configs{Personalized_genome_py}";
        $cmd.= " -i $snpvcf{$sampid}" if($snpvcf{$sampid});
        $cmd.= " -i2 $snvvcf{$sampid}" if($snvvcf{$sampid});
        $cmd.= " -i3 $indelvcf{$sampid}" if($indelvcf{$sampid});
        $cmd.= " -l $configs{pep_length_range}" if($configs{pep_length_range});
        $cmd.= " -r $configs{genome_ref}" if($configs{genome_ref});
        $cmd.= " -s $sampid -o $mydir";
        generateShell($shell, $cmd);
        $scripts{$sampid}{genome} = $shell;
        for my $chr(@Blocks){
            system("mkdir -p $mydir/$chr") unless(-e "$mydir/$chr");
            $shell = "$outdir/$sampid/0.shell/a01.neoantigen/a01.$sampid.$chr.transform.sh";
            $cmd = "if [ -e $mydir/$sampid.$chr.mut.fa.gz ]\nthen\n";
            $cmd.= "\tmv $mydir/$sampid.$chr.mut.fa.gz* $mydir/$chr/ \n";
            $cmd.= "\tmv $mydir/$sampid.$chr.mut.tsv.gz $mydir/$chr/ \n";
            $cmd.= "fi\n";
            $cmd.= "if [ -e $mydir/$chr/$sampid.$chr.mut.tsv.gz ]\nthen\n";
            $cmd.= "\t$configs{orfipy} $mydir/$chr/$sampid.$chr.mut.fa.gz --start $configs{start_condon}";
            $cmd.= " --min $configs{ns_min_length} --outdir $mydir/$chr --dna $sampid.$chr.orf.fa --pep $sampid.$chr.orf_pep.fa";
            $cmd.= " --ignore-case --procs $configs{orfipy_procs_num} && \\\n";
            $cmd.= "\t$configs{orf2mut} $mydir/$chr/$sampid.$chr.orf.fa";
            $cmd.= " $mydir/$chr/$sampid.$chr.mut.tsv.gz $mydir/$chr/$sampid.$chr.filter1.pep.tsv.gz";
            $cmd.= " $mydir/$chr/$sampid.$chr.mut.fa.gz $pep_len[$#pep_len] && \\\n";
            # $cmd.= "\t$configs{PepSplitByLength} $mydir/$chr/$sampid.$chr $mydir/$chr/$sampid.$chr.filter1.pep.tsv.gz\n";
            $cmd.= "\t$configs{PepSplitByLength} $mydir/$chr/$sampid.$chr $BAM{$sampid} $mydir/$chr/$sampid.$chr.filter1.pep.tsv.gz\n";
            $cmd.= "fi";
            generateShell($shell, $cmd);
            $scripts{$sampid}{orf}{$chr} = $shell;
        }
        for my $idx(0..$#pep_len){
            my $len_ns = $ns_len[$idx];
            my $len_pep = $pep_len[$idx];
            $shell = "$outdir/$sampid/0.shell/a01.neoantigen/a01.$sampid.PEP$len_pep.filter.sh";
            # $cmd = "$configs{PepFilterByKmer} $Kmer{$sampid}/$sampid.Tumor.WGS.${len_ns}mer.NT2AA.gz";
            $cmd = "$configs{PepFilterByKmer}";
            $cmd.= " $mydir/$sampid.PEP$len_pep.filter2.tsv $mydir/*/$sampid.*.PEP$len_pep.somatic.tsv";
            generateShell($shell, $cmd);
            $scripts{$sampid}{pep}{$idx} = $shell;
        }
        $shell = "$outdir/$sampid/0.shell/a01.neoantigen/a01.$sampid.Dedup.sh";
        $cmd = "$configs{PepFilterDedup} $CdsMutPro{$sampid}{cds_plus_uniprot_db} $mydir/$sampid.PEP.final.tsv.gz";
        $cmd.= " $mydir/$sampid.uniprot_MutCDS_MutPep.fasta.gz $mydir/$sampid.PEP*.filter2.tsv && \\\n";
        $cmd.= " gzip $mydir/$sampid.PEP*.filter2.tsv && \\\n";
        $cmd.= "gzip -dc $mydir/*/$sampid.*.mut.fa.gz | gzip -c > $mydir/$sampid.mut.fa.gz && \\\n";
        $cmd.= "rm $mydir/*/$sampid.*.mut.fa.gz* && \\\n";
        $cmd.= "gzip -dc $mydir/$Blocks[0]/$sampid.$Blocks[0].mut.tsv.gz | head -1 |gzip -c > $mydir/$sampid.mut.tsv.gz && \\\n";
        $cmd.= "ls $mydir/*/$sampid.*.mut.tsv.gz | while read file; do gzip -dc \$file | tail -n +2 |gzip -c >> $mydir/$sampid.mut.tsv.gz;done && \\\n";
        $cmd.= "cat $mydir/*/$sampid.*.orf.fa | gzip -c > $mydir/$sampid.orf.fa.gz && \\\n";
        $cmd.= "cat $mydir/*/$sampid.*.orf_pep.fa | gzip -c > $mydir/$sampid.orf_pep.fa.gz && \\\n";
        $cmd.= "rm $mydir/*/$sampid.*.mut.tsv.gz\n";
        $cmd.= "rm $mydir/*/$sampid.*.orf.fa $mydir/*/$sampid.*.orf_pep.fa\n";
        $cmd.= "rm $mydir/tmp.*.tsv";
        generateShell($shell, $cmd);
        $scripts{$sampid}{pro_db} = $shell;
        
        $shell = "$outdir/$sampid/0.shell/a01.neoantigen/$sampid.PEP.Annot.sh";
        $cmd = "gzip -dc $mydir/$sampid.mut.tsv.gz |";
        $cmd.= " perl -ne 'next if(/^CHROM/); my \$line = \$_; my (\$CHROM, \$POS, \$REF, \$ALT) = (split /\\t/, \$line, 6)[(0,1,2,4)]; print \"\$CHROM\\t\$POS\\t.\\t\$REF\\t\$ALT\\t.\\tPASS\\t.\\n\"' |";
        $cmd.= " gzip -c > $mydir/$sampid.vcf.gz && \\\n";
        $cmd.= "$configs{ANNOVAR} $mydir/$sampid.vcf.gz $configs{annovar_refdir}";
        $cmd.= " --out $mydir/$sampid.variants.annot --buildver $configs{annovar_buildver} --protocol $configs{annovar_protocol}";
        $cmd.= " --operation $configs{annovar_operation} $configs{annovar_annovar_par} && \\\n";
        $cmd.= "ls $mydir/$sampid.variants.annot* | while read file; do $configs{bgzip} \$file; done && \\\n";
        $cmd.= "perl $configs{PepAnno_pl} $mydir/$sampid.variants.annot.$configs{annovar_buildver}_multianno.vcf.gz $sampid";
        $cmd.= " $mydir/$sampid.PEP.final.tsv.gz $mydir $configs{annovar_protocol}";
        $cmd.= " $configs{gene_alias}" if ($configs{gene_alias});
        $cmd.= " && \\\n";
        $cmd.= "$configs{bgzip} -@ 2 -c $mydir/$sampid.PEP.final.annot.tsv > $mydir/$sampid.PEP.final.tsv.gz\n";
        $cmd.= "rm $mydir/$sampid.vcf.gz $mydir/$sampid.PEP.final.annot.tsv";
        generateShell($shell, $cmd);
        $scripts{$sampid}{pro_anno} = $shell;

    }
}

sub EdgeList{
    my ($day) = @_;
    my $mydir = "$outdir/shell_run";
    system("mkdir -p $outdir/shell_run") unless(-e "$outdir/shell_run");
    open EDGE, ">$outdir/shell_run/edge.$configs{projectname}.$day.list" or die $!;
    for my $sampid(keys %samples){
        print EDGE "$scripts{$sampid}{cds_db}:$configs{cds_db_sge} $scripts{$sampid}{pro_db}:$configs{pro_db_sge}\n";
        for my $chr(@Blocks){
            print EDGE "$scripts{$sampid}{genome}:$configs{genome_sge} $scripts{$sampid}{orf}{$chr}:$configs{orf_sge}\n";
            for my $idx(0..$#pep_len){
                print EDGE "$scripts{$sampid}{orf}{$chr}:$configs{orf_sge} $scripts{$sampid}{pep}{$idx}:$configs{pep_sge}\n";
                print EDGE "$scripts{$sampid}{pep}{$idx}:$configs{pep_sge} $scripts{$sampid}{pro_db}:$configs{pro_db_sge}\n";
            }
        }
        print EDGE "$scripts{$sampid}{pro_db}:$configs{pro_db_sge} $scripts{$sampid}{pro_anno}:$configs{pro_anno_sge}\n";
    }
    close EDGE;
}