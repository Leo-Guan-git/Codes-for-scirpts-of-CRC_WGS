#!/bin/bash
############################################################################################
# Function :MC38 RNA-Seq data preprocessing
# Platform :Linux                                                                          
# Version  :1.0                                                                            
# Author   :YanYX                                                                          
# Date     :2023.03                                                                        
############################################################################################
set -e
workDir="/jdfssz1/ST_SUPERCELLS/P21Z10200N0094/yanyixin/MOUSE/MC38/RNA-Seq"
sampleID="MC38-Pan_merge"
rawdataDir="/jdfssz1/ST_SUPERCELLS/P21Z10200N0094/yanyixin/MOUSE/MC38/RNA-Seq/MC38-Pan_merge/merge_raw"
sampleType="Tumor"
genomeIndexDir="/hwfssz5/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39/STAR_100"

# software
fastpDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/xianghaitao/software/fastp"
rcorrectorDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/xianghaitao/software/Rcorrector/Rcorrector-1.0.5"
STARDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/xianghaitao/software/STAR/STAR-2.7.10a"

# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleID}/${sampleType}/01.DataPreprocess/{fastp,rcorrector}/run

# Step1: Data preprocessing by fastp #
echo "====  Step1 start: Data preprocessing by fastp  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/fastp/run
echo -e "date\n\n${fastpDir}/fastp -i ${rawdataDir}/${sampleID}_1.fq.gz -I ${rawdataDir}/${sampleID}_2.fq.gz -o ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/fastp/${sampleID}.clean.1.fq -O ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/fastp/${sampleID}.clean.2.fq --json ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/fastp/${sampleID}.fastp.json --html ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/fastp/${sampleID}.fastp.html --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --average_qual 10 --length_required 100 --report_title fastp report of ${sampleID}\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > fastp.sh
bash fastp.sh

date
echo "====  Step1 done: Data preprocessing by fastp  ===="

# Step2: Data preprocessing by rcorrector #
echo "====  Step2 start: Data preprocessing by Rcorrector  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/rcorrector/run
echo -e "date\n\nperl ${rcorrectorDir}/run_rcorrector.pl -t 8 -1 ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/fastp/${sampleID}.clean.1.fq   -2 ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/fastp/${sampleID}.clean.2.fq -od  ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/rcorrector\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > rcorrector.sh
bash rcorrector.sh

date
echo "====  Step2 done: Data preprocessing by Rcorrector  ===="

# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleID}/Tumor/02.Aligna2Genome/STAR/run

# Step3: Align to mm30 Genome by STAR #
echo "====  Step3 start: Align to mm30 Genome by STAR  ===="
date

cd ${workDir}/${sampleID}/Tumor/02.Aligna2Genome/STAR/run
echo -e "date\n\n${STARDir}/STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 8 --genomeDir ${genomeIndexDir} --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within KeepPairs --outReadsUnmapped Fastx --readFilesIn ${workDir}/${sampleID}/Tumor/01.DataPreprocess/rcorrector/${sampleID}.clean.1.cor.fq ${workDir}/${sampleID}/Tumor/01.DataPreprocess/rcorrector/${sampleID}.clean.2.cor.fq\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > STAR.sh
bash STAR.sh

date
echo "====  Step3 done: Align to mm30 Genome by STAR  ===="


