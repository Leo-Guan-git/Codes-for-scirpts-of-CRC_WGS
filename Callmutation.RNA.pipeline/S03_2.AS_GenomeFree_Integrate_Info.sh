#!/bin/bash
############################################################################################
# Function :Integrate mutpep & ORF detail info in filter pipeline                          
# Platform :Linux                                                                          
# Version  :1.0                                                                            
# Author   :YanYX                                                                          
# Date     :2023.02                                                                        
############################################################################################
set -e
workDir="/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/yanyixin/MOUSE/MC38/RNA-Seq"
sampleID="MC38_CDX_in-situ_7-Gut"
sampleType="Tumor"
assembleModel="genomeFree"
persegmentpep15Num="1000"

# software
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/bin/bin"
scriptDir="/hwfssz1/ST_SUPERCELLS/P22Z10200N0754/yanyixin/tmp/CRC/denovo/Trinity/Tumor/run/script"

echo "====  Pipeline start:  ===="
date 

# Data preprocessing #
cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter
mkdir -p integrate_info/${assembleModel}/{run,segment_pep}

cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}
for (( i=8;i<=15;i++))
	do
		nohup cut -f 1,2,9,11 ${workDir}/${sampleID}/03.Mutant_analysis_AS/rMATS/${sampleID}.AS.${i}.pep.txt > ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/${sampleID}.AS.${i}.pep.less.txt &
done 

cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/segment_pep
${pythonDir}/python3 ${scriptDir}/file.Segment.8_15pep.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.8.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.9.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.10.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.11.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.12.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.13.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.14.pep ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.15.pep

wait

# Step1: Integrate mutpep & ORF detail info #
echo "====  Step1 start: Integrate mutpep & ORF detail info  ===="
date

## pep15 analysis alone ##
cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/run
for (( i=8;i<=14;i++))
	do
		echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/make.detail.info.AS.v2.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/segment_pep/${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.${i}.pep ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/${sampleID}.AS.${i}.pep.less.txt ${workDir}/${sampleID}/03.Mutant_analysis_AS/rMATS/${sampleID}.AS.GeneInfo.txt ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.pep.result ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/${sampleID}.${sampleType}.${assembleModel}.ORF.detail.less.info ${sampleID}.AS.${sampleType}.pep${i}.Trinity.${assembleModel}.orfipy.detail.info\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > integrate_info.${i}.sh
		nohup bash integrate_info.${i}.sh &
done

ls integrate_info.*.sh | xargs -i echo qsub -cwd -l vf=2G,p=1   -q st.q   -P P21Z10200N0125   {} > run.sh	
	
mkdir pep15

cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/segment_pep
mkdir segment_pep15
cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/segment_pep/segment_pep15

${pythonDir}/python3 ${scriptDir}/file.Segment.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/segment_pep/${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.15.pep ${persegmentpep15Num} ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/segment_pep/segment_pep15 {}/${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.15_{}.pep

ls *.pep | xargs -i echo ${pythonDir}/python3 ${scriptDir}/make.detail.info.AS.v2.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/segment_pep/segment_pep15/{} ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/${sampleID}.AS.${i}.pep.less.txt ${workDir}/${sampleID}/03.Mutant_analysis_AS/rMATS/${sampleID}.AS.GeneInfo.txt ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/match_mutpep/${sampleID}.AS.${sampleType}.pep${i}.${assembleModel}.orfipy.pep.result ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleModel}_Mutpep_Filter/${sampleID}.${sampleType}.${assembleModel}.ORF.detail.less.info {}.detail.info > ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/run/pep15/work.sh
segmentNum=`ls -l *.pep | grep "^-" | wc -l`
cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/run/pep15
for (( var=1;var<=${segmentNum};var++))
	do
		sed -n ''$var'p' work.sh > integrate_info.15.${var}.sh
		nohup bash integrate_info.15.${var}.sh &
done

ls integrate_info.15.*.sh | xargs -i echo qsub -cwd -l vf=2G,p=1   -q st.q   -P P21Z10200N0125   {} > run.sh	

wait

sed -i '1d' *.info
segmentNum=`ls -l *.info | grep "^-" | wc -l`
for (( var=1;var<=${segmentNum};var++))
	do
		cat ${sampleID}.${sampleType}.all.Trinity.orfipy.filtered.15_${var}.pep.detail.info >> ${sampleID}.AS.${sampleType}.pep15.Trinity.${assembleModel}.orfipy.detail.info
done

sed -i '1 i MUTPEPid\tpeptideID\tmutAA\tmutAAlength\tmutSiteNum\tmutSites\tstrand\tgeneNum\tgeneID\tgeneName\tmutTypeNum\tmutTypes\tgenomeType\tORFInfo' ${sampleID}.AS.${sampleType}.pep15.Trinity.${assembleModel}.orfipy.detail.info

nohup mv ${sampleID}.AS.${sampleType}.pep15.Trinity.${assembleModel}.orfipy.detail.info ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel} &
cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel}/run
for (( i=8;i<=14;i++))
	do
		nohup mv ${sampleID}.AS.${sampleType}.pep${i}.Trinity.${assembleModel}.orfipy.detail.info ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleModel} &
done

wait

date
echo "====  Step1 done: Integrate mutpep & ORF detail info  ===="



