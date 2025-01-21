#!/bin/bash
############################################################################################
# Function :RNA-Seq mutant analysis pipeline for gene fusion by STAR-Fusion                
# Platform :Linux                                                                          
# Version  :1.0                                                                           
# Author   :YanYX                                                                          
# Date     :2023.02                                                                        
############################################################################################
set -e
workDir="/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/yanyixin/MOUSE/MC38/RNA-Seq"
sampleID="MC38_CDX_in-situ_6-Gut"
sampleType="Tumor"
genomeIndexDir="/hwfssz1/ST_SUPERCELLS/P22Z10200N0754/yanyixin/tmp/CRC/STAR_Fusion/CTAT-Mutation/CTAT_MUTATION_LIB_SUPPLEMENT/mm39/Mouse_GRCm39_M31_CTAT_lib_Nov092022.plug-n-play/ctat_genome_lib_build_dir"
genomefaDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39"
gtfDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39/gtf"

# software
STARFusionDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/xianghaitao/software/STAR-Fusion-master"
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/zhaoyuntong/softwares/miniconda3/bin"
scriptDir="/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/yanyixin/MOUSE/MC38/RNA-Seq/run/03.Mutant_analysis/STAR_Fusion/STAR-Fusin_tmp/script"

# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleID}/${sampleType}/03.Mutant_analysis_GeneFusion/STAR_Fusion/run
cd ${workDir}/${sampleID}/${sampleType}/03.Mutant_analysis_GeneFusion/STAR_Fusion/run

# Step1: Gene fusion by STAR-Fusion  #
echo "====  Step1 start: Gene fusion by STAR-Fusion  ===="
date

echo -e "date\n\n${STARFusionDir}/STAR-Fusion --genome_lib_dir ${genomeIndexDir} --left_fq ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/rcorrector/${sampleID}.clean.1.cor.fq.gz --right_fq ${workDir}/${sampleID}/${sampleType}/01.DataPreprocess/rcorrector/${sampleID}.clean.2.cor.fq.gz --output_dir ${workDir}/${sampleID}/${sampleType}/03.Mutant_analysis_GeneFusion/STAR_Fusion/run --examine_coding_effect --extract_fusion_reads --FusionInspector inspect --denovo_reconstruct --tmpdir ${workDir}/${sampleID}/${sampleType}/03.Mutant_analysis_GeneFusion/STAR_Fusion/run/tmp --STAR_SortedByCoordinate\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > STARFusion.sh
bash STARFusion.sh

date
echo "====  Step1 done: Gene fusion by STAR-Fusion  ===="

# Step2: Integrate STAR-Fusion analysis result   #
echo "====  Step2 start: Integrate STAR-Fusion analysis result  ===="
date 

cd ${workDir}/${sampleID}/${sampleType}/03.Mutant_analysis_GeneFusion/STAR_Fusion/run
for (( i=8;i<=15;i++))
	do
		echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/S01.getFusionGenesAndPeps.v1.py ${workDir}/${sampleID}/${sampleType}/03.Mutant_analysis_GeneFusion/STAR_Fusion/run/star-fusion.fusion_predictions.abridged.coding_effect.tsv ${i} ${gtfDir}/gencode.vM32.chr_patch_hapl_scaff.annotation.rename.gtf ${genomefaDir}/mm39.fa > ${sampleID}.${sampleType}.${i}.pepInfo1.txt\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > S01.GeneFusion.${i}Info.sh
		nohup bash S01.GeneFusion.${i}Info.sh &
done

wait 

for (( i=8;i<=15;i++))
	do
		nohup awk '$3!~"_" {print}' ${sampleID}.${sampleType}.${i}.pepInfo1.txt > ${workDir}/${sampleID}/${sampleType}/03.Mutant_analysis_GeneFusion/STAR_Fusion/${sampleID}.${sampleType}.${i}.pepInfo.txt &
done

wait

date
echo "====  Step1 done: Gene fusion by STAR-Fusion  ===="