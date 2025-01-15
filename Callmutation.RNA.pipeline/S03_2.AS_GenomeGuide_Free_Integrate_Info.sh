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
sampleID="MC38_CDX_in-situ_6-Gut"
sampleType="Tumor"
assembleGenomeGuide="genomeGuide"
assembleGenomeFree="genomeFree"

# software
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/bin/bin"
scriptDir="/hwfssz1/ST_SUPERCELLS/P22Z10200N0754/yanyixin/tmp/CRC/denovo/Trinity/Tumor/run/script"

# Data preprocessing #
cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info
mkdir -p mergeInfo_genomeGuide_genomeFree/run

# Step1: Integrate genomeGuide & genomeFree detail info #
echo "====  Step1 start: Integrate genomeGuide & genomeFree detail info  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run
for (( i=8;i<=15;i++))
	do
		echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/step1.mergeORF.genomoguide.free.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleGenomeGuide}/${sampleID}.AS.${sampleType}.pep${i}.Trinity.${assembleGenomeGuide}.orfipy.detail.info ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/${assembleGenomeFree}/${sampleID}.AS.${sampleType}.pep${i}.Trinity.${assembleGenomeFree}.orfipy.detail.info ${sampleID}.AS.${sampleType}.pep${i}.Trinity.orfipy.merge.detail.info\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > mergeInfo_S1.${i}.sh
		nohup bash mergeInfo_S1.${i}.sh &
done

wait

for (( i=8;i<=15;i++))
	do
		nohup sed -i "s/nan/ /g" ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.pep${i}.Trinity.orfipy.merge.detail.info &
done

${pythonDir}/python3 ${scriptDir}/generate.pep.ID.AS.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.ID.pep

wait

for (( i=8;i<=15;i++))
	do
		echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/step2.mergeORF.genomoguide.free.AS.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.pep${i}.Trinity.orfipy.merge.detail.info ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.ID.pep ${sampleID}.AS.${sampleType}.pep${i}.Trinity.orfipy.merge.detail.less.info\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > mergeInfo_S2.${i}.sh
		nohup bash mergeInfo_S2.${i}.sh &
done

wait

for (( i=9;i<=15;i++))
	do
		nohup sed -i '1d' ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.pep${i}.Trinity.orfipy.merge.detail.less.info &
done

wait

for (( i=8;i<=15;i++))
	do
		cat ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.pep${i}.Trinity.orfipy.merge.detail.less.info >> ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.unsort.info
done

sort -t '.' -k 3 -n ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.unsort.info > ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.nostartcoden.info

echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/add.unique.startcode.AS.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.nostartcoden.info ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.info\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > add.startcodon.sh
bash add.startcodon.sh

mv ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.info ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree
cp ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.info ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter
nohup rm ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.unsort.info &
nohup rm ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/integrate_info/mergeInfo_genomeGuide_genomeFree/run/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.merger.detail.nostartcoden.info &

wait

echo "====  Step2 done: Integrate mutpep & ORF detail info  ===="
date







