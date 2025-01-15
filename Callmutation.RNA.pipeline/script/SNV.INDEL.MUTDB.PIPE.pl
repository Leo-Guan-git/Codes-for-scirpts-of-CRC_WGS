#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);

=head1 Description

=head1 Version

    Version: 0.1    Date: 2022-12-15   Author: xianghaitao@genomics.cn

=head1 Options

    -vcfSNV         <s> : path to vcf[.gz] file of SNV mutation (mandatory, compatible with gz compressed file)
    -vcfINDEL       <s> : path to vcf[.gz] file of INDEL mutation [optional, compatible with gz compressed file, NOT provided by default]
    -inGenomeFA     <s> : path to reference genome sequece file in fasta format (mandatory, compatible with gz compressed file)
    -inGenomeGTF    <s> : path to reference genome annotation file in gtf format (mandatory, compatible with gz compressed file)
    -proDBfa        <s> : path to reference protein amino acid sequence in fasta format (mandatory, compatible with gz compressed file)
    -ORFdetail      <s> : path to detailed information of predicted ORF from RNAseq data (mandatory)
    -ANNOdir        <s> : path to "ANNOVAR" software directory [optional, default "$Bin/ANNOVAR"]
    -project        <s> : project ID use to submit run job [optional, default "P22Z10200N0754"]
    -mutRegionTag   <s> : mutant region tag use to generate mutate peptide [optional, default "all",
                          assign tags presented in the "Func.refGene" column of "-inAnno" file, multiple region tag can be joined by comma]
    -minPEPlen      <n> : minimum length of mutant peptide [optional, default 8 aa]
    -maxPEPlen      <n> : minimum length of mutant peptide [optional, default 15 aa]
    -splitFileNum   <n> : the number of subfiles the input file will be split into [used to cut file before split run, default: 100]
    -outdir         <s> : path to output directory [optional, default `pwd`]
    -prefix         <s> : prefix of output file [optional, deduce from the basename of '-vcfSNV' by default]
    -help               : show this help information and exit

=cut

my ($vcfSNV, $vcfINDEL, $inGenomeFA, $inGenomeGTF, $proDBfa, $ORFdetail, $ANNOdir, $project, $mutRegionTag, $minPEPlen, $maxPEPlen, $splitFileNum, $outdir, $prefix, $help);
GetOptions (
    "vcfSNV:s" => \$vcfSNV,
    "vcfINDEL:s" => \$vcfINDEL,
    "inGenomeFA:s" => \$inGenomeFA,
    "inGenomeGTF:s" => \$inGenomeGTF,
    "proDBfa:s" => \$proDBfa,
    "ORFdetail:s" => \$ORFdetail,
    "ANNOdir:s" => \$ANNOdir,
    "project:s" => \$project,
    "mutRegionTag:s" => \$mutRegionTag,
    "minPEPlen:n" => \$minPEPlen,
    "maxPEPlen:n" => \$maxPEPlen,
    "splitFileNum:n" => \$splitFileNum,
    "outdir:s" => \$outdir,
    "prefix:s" => \$prefix,
    "help" => \$help
);

die `pod2text $0` if (!defined $vcfSNV || $help);

$ANNOdir ||= "$Bin/ANNOVAR";
$project ||= "P22Z10200N0754";
$mutRegionTag ||= "all";
$minPEPlen ||= 8;
$maxPEPlen ||= 15;
$splitFileNum ||= 100;

$outdir ||= `pwd`;
$outdir = abs_path($outdir);
chomp $outdir;

$prefix ||= basename($vcfSNV);
$prefix =~ s/\.gz//;
$prefix =~ s/\.vcf//;

################### STEP1: annotation mutant sites by ANNOVAR
my $script = "date\n";
$script .= "perl $ANNOdir/table_annovar.pl   $vcfSNV   $ANNOdir/humandb   -buildver hg38   -out ./$prefix.SNV   -protocol refGene   -operation g   -nastring .   -vcfinput   --thread 30   --maxgenethread 30   -polish\n\n";
$script .= "cp ./$prefix.SNV.hg38_multianno.txt  ./$prefix.hg38_multianno.txt\n";
if ($vcfINDEL) {
    $script .= "\nperl $ANNOdir/table_annovar.pl   $vcfINDEL   $ANNOdir/humandb   -buildver hg38   -out ./$prefix.INDEL   -protocol refGene   -operation g   -nastring .   -vcfinput   --thread 30   --maxgenethread 30   -polish\n\n";
    $script .= "sed -n '2,\$'p ./$prefix.INDEL.hg38_multianno.txt >> ./$prefix.hg38_multianno.txt\n";
}
$script .= "perl $Bin/ANNOVARmultiannoMutRegionStat.pl   ./$prefix.hg38_multianno.txt   ./$prefix.hg38_multianno.mutRegion.stat.tsv\n";
$script .= "\nrm -rf ./*.SNV.avinput  ./*.INDEL.avinput  ./*.SNV.hg38_multianno.*  ./*.INDEL.hg38_multianno.*  ./*.refGene.*\n";
$script .= "\ndate\n";

my $id = 1;
my $step = sprintf("%02d", $id);
my $stepSH = "$outdir/step$step\_annoByANNOVAR.sh";
open (SH, "> $stepSH") || die $!;
print SH $script;
close SH;
system("chmod 755 $stepSH");

################### STEP2: perform 6-frame translation around the mutation site to obtain the original mutant peptides
$script = "date\n";
$script .= "perl $Bin/ANNOVARmultianno2MutantPeptides.pl   -inAnno ./$prefix.hg38_multianno.txt   -inGenomeFA $inGenomeFA   -inGenomeGTF $inGenomeGTF   -mutRegionTag $mutRegionTag   -minPEPlen $minPEPlen   -maxPEPlen $maxPEPlen   -outdir ./   -prefix $prefix\n";
$script .= "\ndate\n";

$id ++;
$step = sprintf("%02d", $id);
$stepSH = "$outdir/step$step\_muteSite2mutePep.sh";
open (SH, "> $stepSH") || die $!;
print SH $script;
close SH;
system("chmod 755 $stepSH");

################### STEP3: remove mutant peptides that can be matched in the reference protein database
my $splitDir = "$outdir/splitRun";
system("mkdir -p $splitDir");

$script = "date\n";
$script .= "perl $Bin/file.cutting.pl   -input ./$prefix.SnvIndel.mutant.peptides.raw.gz   -outdir $splitDir   -prefix $prefix.SnvIndel.mutant.peptides.raw   -header yes   -splitFileNum $splitFileNum\n\n";
$script .= "rm -rf $splitDir/qsub.sh\n";
$script .= "for i in \`ls $splitDir/$prefix.SnvIndel.mutant.peptides.raw.*.txt\`\ndo\n";
$script .= "input=\`basename \$i\`\n";
$script .= "fileID=\`basename \$i | awk -F '.' '{print \$(NF-1)}'\`\n";
$script .= "echo \"date\n";
$script .= "perl $Bin/removeMutPepInUniProt.pl   ./\$input   $proDBfa   ./   $prefix.\$fileID   > ./filterMutPepInUniProt.\$fileID.log\n";
$script .= "date\" > $splitDir/$prefix.\$fileID.run.sh\n";
$script .= "chmod 755 $splitDir/$prefix.\$fileID.run.sh\n";
$script .= "echo \"qsub   -cwd   -l vf=1G,num_proc=1   -binding linear:1   -q st.q   -P $project   $prefix.\$fileID.run.sh\" >> $splitDir/qsub.sh\n";
$script .= "done\nchmod 755 $splitDir/qsub.sh\n";

$script .= "echo \"date\n";
$script .= "oneFile=\\`ls ./$prefix.*.SnvIndel.mutant.peptides.noDBmatch.gz | head -1\\`\n";
$script .= "gzip -cd \\\$oneFile | head -1 | gzip > ../$prefix.SnvIndel.mutant.peptides.noDBmatch.gz\n";
$script .= "for i in \\`ls ./$prefix.*.SnvIndel.mutant.peptides.noDBmatch.gz\\`\n";
$script .= "do\n";
$script .= "gzip -cd \\\$i | sed -n '2,\$'p | gzip >> ../$prefix.SnvIndel.mutant.peptides.noDBmatch.gz\n";
$script .= "done\n";
$script .= "\nrm -rf ./$prefix.SnvIndel.mutant.peptides.raw.*.txt  ./$prefix.*.SnvIndel.mutant.peptides.noDBmatch.gz  ./$prefix.*.SnvIndel.noDBmatch.mutPepCount.stat.tsv  ./filterMutPepInUniProt.*.log   ./$prefix.*.run.sh*  ./qsub.sh\n";
$script .= "date\" > $splitDir/merge_removeMutPepInUniProtResults.sh\n";
$script .= "chmod 755 $splitDir/merge_removeMutPepInUniProtResults.sh\n";
$script .= "\ndate\n";

$id ++;
$step = sprintf("%02d", $id);
$stepSH = "$outdir/step$step\_removeMutPepInUniProt.sh";
open (SH, "> $stepSH") || die $!;
print SH $script;
close SH;
system("chmod 755 $stepSH");

################### STEP4: extract mutant peptides that can be matched to any of the ORF predicted from RNAseq data
$script = "date\n";
$script .= "perl $Bin/file.cutting.pl   -input ./$prefix.SnvIndel.mutant.peptides.noDBmatch.gz   -outdir $splitDir   -prefix $prefix.SnvIndel.mutant.peptides.noDBmatch   -header yes   -splitFileNum $splitFileNum\n\n";
$script .= "rm -rf $splitDir/qsub.sh\n\n";
$script .= "for i in \`ls $splitDir/$prefix.SnvIndel.mutant.peptides.noDBmatch.*.txt\`\ndo\n";
$script .= "input=\`basename \$i\`\n";
$script .= "fileID=\`basename \$i | awk -F '.' '{print \$(NF-1)}'\`\n";
$script .= "echo \"date\n";
$script .= "perl $Bin/extractMutPepInORF.pl   ./\$input   $ORFdetail   ./   $prefix.\$fileID\n";
$script .= "date\" > $splitDir/$prefix.\$fileID.run.sh\n\n";
$script .= "chmod 755 $splitDir/$prefix.\$fileID.run.sh\n\n";
$script .= "echo \"qsub   -cwd   -l vf=1G,num_proc=8   -binding linear:8   -q st.q   -P $project   $prefix.\$fileID.run.sh\" >> $splitDir/qsub.sh\n";
$script .= "done\n\nchmod 755 $splitDir/qsub.sh\n\n";

$script .= "echo \"date\n";
$script .= "ls ./$prefix.*.SnvIndel.mutant.peptides.withORFmatch.gz | grep -v part > ./in.list\n\n";
$script .= "perl $Bin/mergeMutPepInORFsplitRunResults.pl   ./in.list   ../   $prefix.SnvIndel.mutant.peptides.noDBmatch\n\n";
$script .= "rm -rf ./in.list  ./$prefix.SnvIndel.mutant.peptides.noDBmatch.*.txt  ./$prefix.*.SnvIndel.mutant.peptides.withORFmatch.gz  ./$prefix.*.run.sh*  ./qsub.sh\n";
$script .= "date\" > $splitDir/merge_extractMutPepInORFresults.sh\n";
$script .= "chmod 755 $splitDir/merge_extractMutPepInORFresults.sh\n\n";

$script .= "echo \"date\n";
$script .= "rm -rf ./runningProgressMonitoring.out.tsv\n";
$script .= "for i in \\`ls -ll ./\*.sh.o\* | sort -k5,5n | awk '{print \\\$NF}'\\`\n";
$script .= "do\nind=\\`basename \\\$i | awk -F '.' '{print \\\$1\\\".\\\"\\\$2}'\\`\n";
$script .= "tag=\\`tail --lines=1 \\\$i | perl -lane '{print join(\\\"\\\\t\\\", \\\$F[1], join(\\\" \\\", \@F[2..\\\$#F]))}'\\`\n";
$script .= "echo -ne \\\"\\\$ind\\\\t\\\$tag\\\\t\\\$i\\\\n\\\" >> ./runningProgressMonitoring.out.tsv\n";
$script .= "done\ndate\" > $splitDir/runningProgressMonitoring.sh\n";
$script .= "chmod 755 $splitDir/runningProgressMonitoring.sh\n";
$script .= "\ndate\n";

$id ++;
$step = sprintf("%02d", $id);
$stepSH = "$outdir/step$step\_extractMutPepInORF.sh";
open (SH, "> $stepSH") || die $!;
print SH $script;
close SH;
system("chmod 755 $stepSH");

################### Remove shorter peptides contained in any of a longer peptides
$script = "date\n";
$script .= "perl $Bin/getPepList4RemoveOverlap.pl   ./$prefix.SnvIndel.mutant.peptides.noDBmatch.withORFmatch.unique.gz   $splitDir   $prefix\n\n";

$script .= "for i in \`ls $splitDir | grep \".vs.\"\`\n";
$script .= "do\nshort=\`echo \$i | awk -F '.' '{print \$1}' | awk -F '_' '{print \$2}'\`\n";
$script .= "perl $Bin/file.cutting.pl   -input $splitDir/\$i/$prefix.length\$short.pep.list   -outdir $splitDir/\$i   -prefix $prefix.length\$short.pep   -header no   -splitFileNum 10\n";
$script .= "done\n\n";

$script .= "rm -rf $splitDir/qsub.sh\n\n";
$script .= "for j in \`ls $splitDir | grep \".vs.\"\`\ndo\n";
$script .= "short=\`echo \$j | awk -F '.' '{print \$1}' | awk -F '_' '{print \$2}'\`\n";
$script .= "long=\`echo \$j | awk -F '.vs.' '{print \$NF}'\`\n";
$script .= "for k in \`ls $splitDir/\$j/$prefix.length\$short.pep.*.txt\`\n  do\n";
$script .= "input=\`basename \$k\`\n";
$script .= "file=\`echo \$input | awk -F '.' '{print \$(NF-1)}'\`\n";
$script .= "echo \"date\nperl $Bin/removeOverlapMutPEP.pl   ./\$input   ./$prefix.length\$long.pep.list   ./$prefix.\$j.\$file.retained.peptides.list  > ./$prefix.\$j.\$file.log\n";
$script .= "date\" > $splitDir/\$j/run.\$j.\$file.sh\n";
$script .= "chmod 755 $splitDir/\$j/run.\$j.\$file.sh\n\n";
$script .= "echo \"qsub   -wd $splitDir/\$j   -l vf=1G,num_proc=1   -binding linear:1   -q st.q   -P $project   $splitDir/\$j/run.\$j.\$file.sh\" >> $splitDir/qsub.sh\n";
$script .= "done\ndone\nchmod 755 $splitDir/qsub.sh\n\n";

$script .= "echo \"date\nls ./length_\*/$prefix.\*.retained.peptides.list > ./retained.peptides.list\n\n";
$script .= "perl $Bin/mergeRemoveOverlapMutPEPresults.pl   ./retained.peptides.list   $outdir/$prefix.SnvIndel.mutant.peptides.noDBmatch.withORFmatch.unique.gz   $outdir/$prefix.SnvIndel.mutant.peptides.noDBmatch.withORFmatch.noOverlap.gz\n\n";
$script .= "rm -rf ./length_*  ./retained.peptides.list\n";
$script .= "date\"> $splitDir/merge_removeOverlapMutPEPresults.sh\n";
$script .= "chmod 755 $splitDir/merge_removeOverlapMutPEPresults.sh\n\n";
$script .= "date\n";

$id ++;
$step = sprintf("%02d", $id);
$stepSH = "$outdir/step$step\_removeOverlapMutPEP.sh";
open (SH, "> $stepSH") || die $!;
print SH $script;
close SH;
system("chmod 755 $stepSH");

################### Generate fasta file for mutant database construction from candidate mutant peptdies featue file
$script = "date\n";
$script .= "perl $Bin/getSNVindelMUTdb.pl   ./$prefix.SnvIndel.mutant.peptides.noDBmatch.withORFmatch.noOverlap.gz   ./   $prefix\n\n";
$script .= "ls ./$prefix.SnvIndel.MUTPEP.fasta   ../alternativeSplicing/$prefix.AS.Tumor.all.Trinity.orfipy.filtered.pep.fa   ../Funsion/$prefix.Fusion.Tumor.all.Trinity.orfipy.filtered.pep.fa > ./$prefix.allType.MUTPEP.list\n";
$script .= "perl $Bin/../merge.SnvIndel.AS.Fusion.MUTPEPdb.pl   ./$prefix.allType.MUTPEP.list   ./$prefix.allType.MUTPEPdb.fasta\n";
$script .= "perl $Bin/MUTPEPfa2MUTPROfa.pl   ./$prefix.allType.MUTPEPdb.fasta   ./$prefix.UniProtWithAllTypeMUTPRO.fasta   ./$prefix.allType.PRO2PEP.position.tsv   $proDBfa\n\n";
$script .= "date\n";

$id ++;
$step = sprintf("%02d", $id);
$stepSH = "$outdir/step$step\_getSNVindelMUTdb.sh";
open (SH, "> $stepSH") || die $!;
print SH $script;
close SH;
system("chmod 755 $stepSH");
