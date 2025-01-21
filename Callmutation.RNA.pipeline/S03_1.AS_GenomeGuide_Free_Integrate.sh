#!/bin/bash
############################################################################################
# Function :Integrate Trinity Genomeguide & Genomefree filtered mutpeps                    
# Platform :Linux                                                                          
# Version  :1.0                                                                            
# Author   :YanYX                                                                          
# Date     :2023.02                                                                        
############################################################################################
set -e
workDir="/jdfssz1/ST_SUPERCELLS/P21Z10200N0094/yanyixin/MOUSE/MC38/RNA-Seq"
sampleID="MC38_CDX_in-situ_5-Gut"
sampleType="Tumor"
assembleGenomeGuide="genomeGuide"
assembleGenomeFree="genomeFree"

# software
pythonDir="/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/yanyixin/bin/bin"
scriptDir="/hwfssz1/ST_SUPERCELLS/P22Z10200N0754/yanyixin/tmp/CRC/denovo/Trinity/Tumor/run/script"


# Data preprocessing #
cd ${workDir}
mkdir -p ${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/{filter_contain,generate_fa}/run

# Step1: Filter by Genomeguide & Genomefree contain #
echo "====  Step1 start: Filter by Genomeguide & Genomefree contain  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/filter_contain/run
echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/filter_contain.py ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleGenomeFree}_Mutpep_Filter/filter_contain/${sampleID}.AS.${sampleType}.all.${assembleGenomeFree}.orfipy.filtered.notcontain.pep ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleGenomeGuide}_Mutpep_Filter/filter_contain/${sampleID}.AS.${sampleType}.all.${assembleGenomeGuide}.orfipy.filtered.notcontain.pep ${sampleID}.AS.${sampleType}.all.${assembleGenomeFree}.notcontain.${assembleGenomeGuide}.orfipy.pep\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > filter_contain.guide_free.sh
bash filter_contain.guide_free.sh

mv ${sampleID}.AS.${sampleType}.all.${assembleGenomeFree}.notcontain.${assembleGenomeGuide}.orfipy.pep ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/filter_contain
cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/filter_contain
cat ${workDir}/${sampleID}/${sampleType}/05.AS_${assembleGenomeGuide}_Mutpep_Filter/filter_contain/${sampleID}.AS.${sampleType}.all.${assembleGenomeGuide}.orfipy.filtered.notcontain.pep ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/filter_contain/${sampleID}.AS.${sampleType}.all.${assembleGenomeFree}.notcontain.${assembleGenomeGuide}.orfipy.pep > ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep
cp ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter

date
echo "====  Step1 done: Filter by Genomeguide & Genomefree contain  ===="

# Step2: Generate mutpep database in fasta format #
echo "====  Step2 start: Generate mutpep database in fasta format  ===="
date

cd ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/generate_fa/run
echo -e "date\n\n${pythonDir}/python3 ${scriptDir}/makefasta.AS.py ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep.fa\n\ndate\n\nqstat -j \$JOB_ID | grep cpu" > generate_fa.sh
bash generate_fa.sh

cp ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep.fa ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter
mv ${sampleID}.AS.${sampleType}.all.Trinity.orfipy.filtered.pep.fa ${workDir}/${sampleID}/${sampleType}/05.AS_Intergrate_genomeGuide_genomeFree_Mutpep_Filter/generate_fa

date
echo "====  Step2 done: Generate mutpep database in fasta format  ===="

