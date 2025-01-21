#!/bin/bash
############################################################################################
# Function :RNA-Seq mutant peptides filter pipeline for Trinity Genomeguide or Genomefree  
# Platform :Linux                                                                          
# Version  :1.1                                                                            
# Author   :YanYX                                                                          
# Date     :2023.02                                                                        
############################################################################################
set -e
workDir="/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/yanyixin/MOUSE/MC38/RNA-Seq"
sampleID="MC38_CDX_in-situ_6-Gut"
sampleType="Tumor"
assembleModel="genomeGuide"
persegmentORFNum="500000"
uniprotdbDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/CRC_merged"

# software
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/bin/bin"
scriptDir="/hwfssz1/ST_SUPERCELLS/P22Z10200N0754/yanyixin/tmp/CRC/denovo/Trinity/Tumor/run/script"


# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter

cut -f 1-9 ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/${sampleID}.${sampleType}.${assembleModel}.ORF.detail.info.txt > ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/${sampleID}.${sampleType}.${assembleModel}.ORF.detail.less.info

# zcat ${rawORFDir}/${sampleID}.Tumor.${assembleModel}.ORF.detail.info.txt.gz | cut -f 1-9 > ${ORFDir}/${assembleModel}/${sampleID}/${sampleID}.Tumor.${assembleModel}.ORF.detail.less.info

cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter
mkdir -p match_mutpep/{run,segment_ORF}
cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/run

for (( i=8;i<=15;i++))
	do
		mkdir -p pep${i}
done

for (( i=8;i<=15;i++))
	do
		nohup awk -F '\t' '{print $2}' ${workDir}/${sampleID}/03.Mutant_analysis_AS/rMATS/${sampleID}.AS.${i}.pep.txt | awk ' !x[$0]++' | sed '/_/d' > ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/${sampleID}.AS.peps.${i}.unique &
done

${pythonDir}/python3 ${scriptDir}/file.Segment.py ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/${sampleID}.${sampleType}.${assembleModel}.orfipy.predicted.ORF.pep.line.fasta ${persegmentORFNum} ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/segment_ORF {}/${sampleID}.${sampleType}.${assembleModel}.orfipy.predicted.ORF.pep.line.{}.fasta

wait

# Step1: Filter by mutpeps match ORFs #
echo "====  Step1 start: Filter by mutpeps match ORFs  ===="
date

for (( i=8;i<=15;i++))
	do 
		cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/segment_ORF
		ls *.fasta | xargs -i echo ${pythonDir}/python3 ${scriptDir}/match.denovo.AA.mutpep.v2.py ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/${sampleID}.AS.peps.${i}.unique ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/segment_ORF/{} {}.match.pep${i}.result > ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/run/pep${i}/work.sh
		segmentNum=`ls -l *.fasta | grep "^-" | wc -l`
		cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/run/pep${i}
		for (( var=1;var<=${segmentNum};var++))
			do
				sed -n ''$var'p' work.sh > match.${i}.${var}.sh
				nohup bash match.${i}.${var}.sh &
		done
done

for (( i=8;i<=15;i++))
	do
		cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/run/pep${i}
		ls match* | xargs -i echo qsub -cwd -l vf=2G,p=1   -q st.q   -P P21Z10200N0125   {} > run.sh
done

wait

for (( i=8;i<=15;i++))
	do 
		cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/run/pep${i}
		cat *.fasta.match.pep${i}.result > ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.pep.result
		nohup awk -F '\t' '{print $1}' ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.pep.result | awk ' !x[$0]++' > ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.unique.pep.result &
done

wait

for (( i=8;i<=15;i++))
	do
		cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/run/pep${i}
		nohup mv ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.pep.result ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep &
done

wait

date
echo "====  Step1 done: Filter by mutpeps match ORFs  ===="

# Step2: Filter by unalign UniProtDB #
echo "====  Step2 start: Filter by unalign UniProtDB  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter
mkdir -p unalign_UniProt/run

cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/unalign_UniProt/run
for (( i=8;i<=15;i++))
	do 
		echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/filter_contain.py ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.unique.pep.result ${uniprotdbDir}/uniprot-filtered-organism_Homo.sapiens_9606_AND_review.line.unique.fasta ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.UniProtfilter.pep\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > unalign.${i}.sh
		nohup bash unalign.${i}.sh &
done

wait

nohup cat ${sampleID}.AS.${sampleType}.pep9.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep10.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep11.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep12.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep13.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep > ${sampleID}.AS.${sampleType}.pep9_15.${assembleModel}.orfipy.UniProtfilter.pep &
nohup cat ${sampleID}.AS.${sampleType}.pep10.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep11.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep12.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep13.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep > ${sampleID}.AS.${sampleType}.pep10_15.${assembleModel}.orfipy.UniProtfilter.pep &
nohup cat ${sampleID}.AS.${sampleType}.pep11.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep12.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep13.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep > ${sampleID}.AS.${sampleType}.pep11_15.${assembleModel}.orfipy.UniProtfilter.pep &
nohup cat ${sampleID}.AS.${sampleType}.pep12.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep13.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep > ${sampleID}.AS.${sampleType}.pep12_15.${assembleModel}.orfipy.UniProtfilter.pep &
nohup cat ${sampleID}.AS.${sampleType}.pep13.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep > ${sampleID}.AS.${sampleType}.pep13_15.${assembleModel}.orfipy.UniProtfilter.pep &
nohup cat ${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep > ${sampleID}.AS.${sampleType}.pep14_15.${assembleModel}.orfipy.UniProtfilter.pep &

wait

for (( i=8;i<=15;i++))
	do
		nohup mv ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.UniProtfilter.pep ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/unalign_UniProt &
done

wait

date
echo "====  Step2 done: Filter by unalign UniProtDB  ===="

# Step3: Filter by 8-15 contain #
echo "====  Step3 start: Filter by 8-15 contain  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter
mkdir -p filter_contain/run

## pep14 analysis alone ##
cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/filter_contain/run
for (( i=8;i<=13;i++))
	do
		echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/filter_contain.py ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/unalign_UniProt/${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.UniProtfilter.pep ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/unalign_UniProt/run/${sampleID}.AS.${sampleType}.pep$(($i+1))_15.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.UniProtfilter.notcontain.pep\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > f_contain.${i}.sh
		nohup bash f_contain.${i}.sh &
done
		
echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/filter_contain.py ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/unalign_UniProt/${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.pep ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/unalign_UniProt/${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep ${sampleID}.AS.${sampleType}.pep14.${assembleModel}.orfipy.UniProtfilter.notcontain.pep\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > f_contain.14.sh
nohup bash f_contain.14.sh &

wait

for (( i=8;i<=14;i++))
	do
		nohup mv ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.UniProtfilter.notcontain.pep ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/filter_contain &
done

cp ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/unalign_UniProt/${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.pep ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/filter_contain/${sampleID}.AS.${sampleType}.pep15.${assembleModel}.orfipy.UniProtfilter.notcontain.pep

wait

cd ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/filter_contain
for (( i=8;i<=15;i++))
	do
		cat ${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.UniProtfilter.notcontain.pep >> ${sampleID}.AS.${sampleType}.all.${assembleModel}.orfipy.filtered.notcontain.pep
done

date
echo "====  Step3 done: Filter by 8-15 contain  ===="