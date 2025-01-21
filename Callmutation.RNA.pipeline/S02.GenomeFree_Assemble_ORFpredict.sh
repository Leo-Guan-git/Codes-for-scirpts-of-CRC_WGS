#!/bin/bash
###################################################################################################
# Function :RNA-Seq data assemble by Trinity and ORF predict by orfipy                            
# Platform :Linux                                                                                 
# Version  :1.0                                                                                   
# Author   :YanYX                                                                                 
# Date     :2023.03                                                                               
###################################################################################################
set -e
workDir="/jdfssz1/ST_SUPERCELLS/P21Z10200N0094/yanyixin/MOUSE/MC38/RNA-Seq"
sampleID="MC38_CDX_in-situ_2-Gut_merge"
sampleType="Tumor"
assembleModel="genomeFree"

# software
perlDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/zhaoyuntong/softwares/miniconda3/bin"
TrinityDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/xianghaitao/software/Trinity/trinityrnaseq-v2.14.0"
orfipyDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/zhaoyuntong/softwares/miniconda3/bin"

# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/run

# Step1: Assemble by Trinity #
echo "====  Step1 start: Assemble by Trinity  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/run
echo -e "date\n\n${perlDir}/perl /hwfssz1/ST_SUPERCELLS/P22Z10200N0754/xianghaitao/NEOANTIGEN/bin/in_development/changeSTARunalignedFQformat.pl ${workDir}/${sampleID}/${sampleType}/02.Aligna2Genome/STAR/run/Unmapped.out.mate1 ${workDir}/${sampleID}/${sampleType}/02.Aligna2Genome/STAR/run/Unmapped.out.mate2 ${workDir}/${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/run ${sampleID}.${sampleType}.${assembleModel}\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > changeSTARunalignedFQformat.sh
bash changeSTARunalignedFQformat.sh

echo -e "date\n\n${TrinityDir}/Trinity --seqType fq --left ${workDir}/${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/run/${sampleID}.${sampleType}.${assembleModel}.1.fq --right ${workDir}/${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/run/${sampleID}.${sampleType}.${assembleModel}.2.fq --CPU 8 --max_memory 50G --output ${workDir}/${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/trinity_out\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > ${assembleModel}_Assemby_Trinity.sh
bash ${assembleModel}_Assemby_Trinity.sh

date
echo "====  Step1 done: Assemble by Trinity  ===="

# Step2: ORF predict by orfipy #
echo "====  Step2 start: ORF predict by orfipy  ===="
date

cd ${workDir}
mkdir -p ${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/run
cd ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/run

echo -e "date\n\n${orfipyDir}/orfipy ${workDir}/${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/trinity_out.Trinity.fasta --start ATG,ATA,ATT,ATC,AAG,ACG,AGG,TTG,CTG,GTG --stop TAA,TAG,TGA --min 24 --outdir ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy --dna "${sampleID}.${sampleType}.${assembleModel}.orfipy.predicted.ORF.NT.fasta"\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > ${assembleModel}_ORFpredict_orfipy.sh
nohup bash ${assembleModel}_ORFpredict_orfipy.sh &

echo -e "date\n\n${orfipyDir}/orfipy ${workDir}/${sampleID}/${sampleType}/03.${assembleModel}_Assemby/Trinity/trinity_out.Trinity.fasta --start ATG,ATA,ATT,ATC,AAG,ACG,AGG,TTG,CTG,GTG --stop TAA,TAG,TGA --min 24 --outdir ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy --pep "${sampleID}.${sampleType}.${assembleModel}.orfipy.predicted.ORF.pep.fasta"\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > ${assembleModel}_ORFpredict_pep_orfipy.sh
bash ${assembleModel}_ORFpredict_pep_orfipy.sh

wait

echo -e "date\n\n${perlDir}/perl /hwfssz1/ST_SUPERCELLS/P22Z10200N0754/xianghaitao/NEOANTIGEN/bin/in_development/predictedORFsummary.pl ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/${sampleID}.${sampleType}.${assembleModel}.orfipy.predicted.ORF.NT.fasta ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy ${sampleID}.${sampleType}.${assembleModel}\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > ${assembleModel}_ORFpredict_detail_info.sh
bash ${assembleModel}_ORFpredict_detail_info.sh

# gzip ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/${sampleID}.${sampleType}.${assembleModel}.ORF.detail.info.txt

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/${sampleID}.${sampleType}.${assembleModel}.orfipy.predicted.ORF.pep.fasta > ${workDir}/${sampleID}/${sampleType}/04.${assembleModel}_ORF_Predict/orfipy/${sampleID}.${sampleType}.${assembleModel}.orfipy.predicted.ORF.pep.line.fasta

wait

date
echo "====  Step2 done: ORF predict by orfipy  ===="





