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
sampleIDNormal="MC38_CDX_in-situ_10-Gut"
sampleIDTumor="MC38_CDX_in-situ_10-Gut"
sampleTypeNormal="Normal"
sampleTypeTumor="Tumor"

# software
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/bin/bin"
scriptDir="/hwfssz1/ST_SUPERCELLS/P22Z10200N0754/yanyixin/tmp/CRC/denovo/Trinity/Tumor/run/script"

# Step3: Integrate RSEM quantitative result #
echo "====  Step3 start: Integrate RSEM quantitative result  ===="
date

cd ${workDir}/${sampleIDTumor}
mkdir -p 03.Quantitative_Intergrate_Tumor_Normal/run
cd ${workDir}/${sampleIDTumor}/03.Quantitative_Intergrate_Tumor_Normal/run

echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/step2.quantitative.result.integrate.v2.py --input1 ${workDir}/${sampleIDTumor}/${sampleTypeTumor}/03.Quantitative_analysis/RSEM/${sampleIDTumor}.${sampleTypeTumor}.quantitative.integrate.txt --input2 ${workDir}/${sampleIDNormal}/${sampleTypeNormal}/03.Quantitative_analysis/RSEM/${sampleIDNormal}.${sampleTypeNormal}.quantitative.integrate.txt --prefix ${sampleIDTumor} --outdir ${workDir}/${sampleIDTumor}/03.Quantitative_Intergrate_Tumor_Normal\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > step2.quantitative.result.integrate.sh
bash step2.quantitative.result.integrate.sh

date
echo "====  Step3 done: Integrate RSEM quantitative result  ===="
