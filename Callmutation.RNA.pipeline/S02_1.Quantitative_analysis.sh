#!/bin/bash
###################################################################################################
# Function :RNA-Seq data quantitative analysis                                                    
# Platform :Linux                                                                                 
# Version  :1.0                                                                                   
# Author   :YanYX                                                                                 
# Date     :2023.07                                                                               
###################################################################################################
set -e
workDir="/hwfssz5/ST_SUPERCELLS/P21Z10200N0125/yanyixin/Mouse/RNA-Seq"
sampleID="MC38_CDX_in-situ_10-Gut"
sampleType="Tumor"
RSEMIndexDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/CRC_merged/genome_index/RSEM/mm39"
RSEMIndexFile="mm39"
gftDir="/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/yanyixin/MOUSE/MC38/RNA-Seq"
gtfFile="gencode.vM32.chr_patch_hapl_scaff.annotation.rename.less.gtf"
# RSEMIndexDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/CRC_merged/genome_index/RSEM/hg38"
# RSEMIndexFile="hg38"
# gftDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/CRC_merged"
# gtfFile="gencode.v38.chr_patch_hapl_scaff.annotation.rename.less.gtf"

# software
RSEMDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/software/RSEM/RSEM-1.3.3"
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/bin/bin"
scriptDir="/hwfssz1/ST_SUPERCELLS/P22Z10200N0754/yanyixin/tmp/CRC/denovo/Trinity/Tumor/run/script"

# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleID}/${sampleType}/03.Quantitative_analysis/RSEM/run

# Step1: Quantitative analysis by RSEM #
echo "====  Step1 start: Quantitative analysis by RSEM  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/03.Quantitative_analysis/RSEM/run
echo -e "date\n\n${RSEMDir}/rsem-calculate-expression --paired-end --alignments ${workDir}/${sampleID}/${sampleType}/02.Aligna2Genome/STAR/run/Aligned.toTranscriptome.out.bam ${RSEMIndexDir}/${RSEMIndexFile} ${workDir}/${sampleID}/${sampleType}/03.Quantitative_analysis/RSEM/${sampleID}\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > RSEM.calculate.expression.sh
bash RSEM.calculate.expression.sh

date
echo "====  Step1 done: Quantitative analysis by RSEM  ===="

# Step2: Integrate RSEM quantitative result #
echo "====  Step2 start: Integrate RSEM quantitative result  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/03.Quantitative_analysis/RSEM/run
echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/step1.quantitative.result.integrate.v2.py --input ${workDir}/${sampleID}/${sampleType}/03.Quantitative_analysis/RSEM/${sampleID}.genes.results --gtf ${gftDir}/${gtfFile} --prefix ${sampleID} --type ${sampleType} --outdir ${workDir}/${sampleID}/${sampleType}/03.Quantitative_analysis/RSEM\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > step1.quantitative.result.integrate.sh
bash step1.quantitative.result.integrate.sh

date
echo "====  Step2 done: Integrate RSEM quantitative result  ===="
