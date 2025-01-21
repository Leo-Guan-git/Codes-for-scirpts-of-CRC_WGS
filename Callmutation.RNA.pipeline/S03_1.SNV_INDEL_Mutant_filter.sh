#!/bin/bash
###################################################################################################
# Function :RNA-Seq data SNV/InDel mutant peptide filter pipeline                                 
# Platform :Linux                                                                                 
# Version  :1.0                                                                                   
# Author   :YanYX                                                                                 
# Date     :2023.03                                                                               
###################################################################################################

dir="/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/yanyixin/MOUSE/MC38/RNA-Seq"
sample="MC38_CDX_in-situ_7-Gut"
sampleType="Tumor"
assembleGenomeGuide="genomeGuide"
assembleGenomeFree="genomeFree"
SNVprefix="${sample}_SNV"
INDELprefix="${sample}_INDEL"
vcfDir=/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/yanyixin/MOUSE/MC38/RNA-Seq/VCF
genomeFA=/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39/mm39.fa
genomeGTF=/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39/gtf/gencode.vM32.chr_patch_hapl_scaff.annotation.rename.gtf
refAA=/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39/uniprot-UP000000589_10090_includeIsoform_tr-2022.11.30.fasta.gz

# software
perlDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/zhaoyuntong/softwares/miniconda3/bin"


# Data preprocessing #
mkdir -p ${dir}/predictedORF
sed -n '1!p' ${dir}/${sample}/${sampleType}/04.${assembleGenomeGuide}_ORF_Predict/orfipy/${sample}.${sampleType}.${assembleGenomeGuide}.ORF.detail.info.txt > ${dir}/predictedORF/${sample}.${sampleType}.${assembleGenomeGuide}.ORF.detail.info.txt

cat ${dir}/${sample}/${sampleType}/04.${assembleGenomeFree}_ORF_Predict/orfipy/${sample}.${sampleType}.${assembleGenomeFree}.ORF.detail.info.txt ${dir}/predictedORF/${sample}.${sampleType}.${assembleGenomeGuide}.ORF.detail.info.txt > ${dir}/predictedORF/${sample}.ORF.detail.info

gzip ${dir}/predictedORF/${sample}.ORF.detail.info

rm ${dir}/predictedORF/${sample}.${sampleType}.${assembleGenomeGuide}.ORF.detail.info.txt

# Step1: Mutant peptide Filter #
echo "====  Step1 start: Mutant peptide Filter  ===="
date

cd ${dir}/${sample}/${sampleType}
mkdir -p 05.SNV_INDEL_Mutpep_Filter/run
cd ${dir}/${sample}/${sampleType}/05.SNV_INDEL_Mutpep_Filter/run

echo -e "date\n\n${perlDir}/perl /hwfssz1/ST_SUPERCELLS/P22Z10200N0754/xianghaitao/NEOANTIGEN/bin/mutationPeptideDBbuild/SNV_INDEL/SNV.INDEL.MUTDB.PIPE.pl -vcfSNV $vcfDir/SNV/$sample.merged_snv.vcf.gz -vcfINDEL $vcfDir/INDEL/$sample.merged_indel.vcf.gz -inGenomeFA $genomeFA -inGenomeGTF $genomeGTF -proDBfa $refAA -ORFdetail ${dir}/predictedORF/$sample.ORF.detail.info.gz -splitFileNum 150 -outdir ${dir}/${sample}/${sampleType}/05.SNV_INDEL_Mutpep_Filter/run -prefix $sample\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > SNV_INDEL_Filter_pipeline.sh
bash SNV_INDEL_Filter_pipeline.sh

date
echo "====  Step1 done: Mutant peptide Filter  ===="
