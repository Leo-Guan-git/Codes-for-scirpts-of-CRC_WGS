#!/bin/bash
############################################################################################
# Function :RNA-Seq mutant analysis pipeline for alternative splicing by rMATS             
# Platform :Linux                                                                          
# Version  :1.0                                                                            
# Author   :YanYX                                                                          
# Date     :2023.02                                                                        
############################################################################################
set -e
workDir="/jdfssz1/ST_SUPERCELLS/P21Z10200N0094/yanyixin/MOUSE/MC38/RNA-Seq"
sampleIDNormal="MC38_CDX_in-situ_2-Gut_merge"
sampleIDTumor="MC38_CDX_in-situ_2-Gut_merge"
sampleTypeb1="Normal"
sampleTypeb2="Tumor"
genomefaDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39"
genomefaFile="mm39.fa"
gtfDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/01.database/mm39/gtf"
gtfFile="gencode.vM32.chr_patch_hapl_scaff.annotation.rename.gtf"

# software
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/zhaoyuntong/softwares/miniconda3/bin"
rMATSDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/zhaoyuntong/softwares/miniconda3/bin"


# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run
cd ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run
ls ${workDir}/${sampleIDNormal}/${sampleTypeb1}/02.Aligna2Genome/STAR/run/Aligned.sortedByCoord.out.bam > ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run/b1.list
ls ${workDir}/${sampleIDTumor}/${sampleTypeb2}/02.Aligna2Genome/STAR/run/Aligned.sortedByCoord.out.bam > ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run/b2.list

# Step1: Alternative splicing by rMATS #
echo "====  Step1 start: Alternative splicing by rMATS  ===="
date

echo -e "date\n\n${pythonDir}/python3 ${rMATSDir}/rmats.py --b1 ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run/b1.list --b2 ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run/b2.list --gtf ${gtfDir}/${gtfFile} -t paired --readLength 100 --nthread 12 --od ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run --tmp ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run/rmats\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > rMATS.sh
bash rMATS.sh

date
echo "====  Step1 done: Alternative splicing by rMATS  ===="

# Step2: Integrate rMATS analysis result #
echo "====  Step2 start: Integrate rMATS analysis result   ===="
date 

cd ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run
for ASType in {SE,RI,MXE,A5SS,A3SS}
	do
		echo -e "date\n\n${pythonDir}/python3 /jdfssz1/ST_SUPERCELLS/P21Z10200N0094/yanyixin/MOUSE/MC38/RNA-Seq/run/03.Mutant_analysis/AS_rMATS/rmats_tmp/script/S01.parseAS.v1.py ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run/${ASType}.MATS.JC.txt ${ASType} ${genomefaDir}/${genomefaFile} > ${ASType}.MATS.JC.Info.txt\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > S01.AS.${ASType}Info.sh
		nohup bash S01.AS.${ASType}Info.sh &
done

wait

for ASType in {A3SS,A5SS,MXE,RI,SE}
	do
		sed -i 's/"//g' ${ASType}.MATS.JC.Info.txt
		cat ${ASType}.MATS.JC.Info.txt >> ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/${sampleIDTumor}.AS.GeneInfo.txt
done

for ASType in {SE,RI,MXE,A5SS,A3SS}
	do	
		for (( i=8;i<=15;i++))
			do
				echo -e "date\n\n${pythonDir}/python3 /jdfssz1/ST_SUPERCELLS/P21Z10200N0094/yanyixin/MOUSE/MC38/RNA-Seq/run/03.Mutant_analysis/AS_rMATS/rmats_tmp/script/S02.steping_trans.v2.py ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/run/${ASType}.MATS.JC.Info.txt ${i} ${ASType}.MATS.JC.${i}.Info.txt\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > S02.AS.${ASType}.${i}Info.sh
				echo qsub -cwd -l vf=2G,p=1   -q st.q   -P P21Z10200N0125   S02.AS.${ASType}.${i}Info.sh >> run.sh
				bash S02.AS.${ASType}.${i}Info.sh
		done
done

for (( i=8;i<=15;i++))
	do	
		for ASType in {A3SS,A5SS,MXE,RI,SE}
			do
				cat  ${ASType}.MATS.JC.${i}.Info.txt >> ${workDir}/${sampleIDTumor}/03.Mutant_analysis_AS/rMATS/${sampleIDTumor}.AS.${i}.pep.txt
		done
done

date
echo "====  Step2 done: Integrate rMATS analysis result  ===="