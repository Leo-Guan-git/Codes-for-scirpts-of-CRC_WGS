import os
import pandas as pd
import datetime

""" Function """
# creare pepline work dir [result,raw]
def create_pepline_workDir(infile, outdir):
    sampleInfo_data = pd.read_csv(infile, sep = '\t', header=None, names=['Sample', 'Type', 'anno', 'psmFile'])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    sampleInfo_data['ID'] = sampleInfo_data['Sample'] + '#' +  sampleInfo_data['Type'] + '#' + sampleInfo_data['anno']
    sampleInfo_dict = sampleInfo_data.set_index("ID").to_dict(orient="index") # HCT116MG_D181211#tumor#SEC  wildcards.sample
    # create sample dir[easy for qsub]
    for sample in sampleInfo_dict.keys():
        sampleDir = os.path.join(outdir, sample)
        if not os.path.exists(sampleDir):
            os.makedirs(sampleDir) 
    return sampleInfo_dict


""" main """
infile = config['selfSet']['infile']
outdir = os.path.abspath(config["selfSet"]["outdir"])
sampleInfo_dict = create_pepline_workDir(infile, outdir)




rule all:
    input:
        ## spectra similarity(SpecSim)
        raw_format_mgf = expand(os.path.join(outdir, "{sample}", "{sample}.raw.format.mgf"), sample=sampleInfo_dict.keys()),
        pDeep2_input = expand(os.path.join(outdir, "{sample}", "{sample}.pDeep2_input.txt"), sample=sampleInfo_dict.keys()),
        pDeep2_output = expand(os.path.join(outdir, "{sample}", "{sample}.pDeep2_output.mgf"), sample=sampleInfo_dict.keys()),
        pDeep2_format_mgf = expand(os.path.join(outdir, "{sample}", "{sample}.pDeep2_output.format.mgf"), sample=sampleInfo_dict.keys()),
        similarity_result = expand(os.path.join(outdir, "{sample}", "{sample}.rawVSpDeep.spectra.similarity.txt"), sample=sampleInfo_dict.keys()),
        sample_psm_result = expand(os.path.join(outdir, "{sample}", "{sample}.addSpectraSimilarityScore.PSM.txt"), sample=sampleInfo_dict.keys()),
        SpecSim_aggregate_file = os.path.join(outdir, config["selfSet"]["prefix"] + ".addSpectraSimilarityScore.PSM.txt"),
        ## retetion time
        AutoRT_testData = expand(os.path.join(outdir, "{sample}", "{sample}.AutoRTtestData.txt"), sample=sampleInfo_dict.keys()),
        DeepLC_testData = expand(os.path.join(outdir, "{sample}", "{sample}.DeepLCpredData.csv"), sample=sampleInfo_dict.keys()),
        AutoRT_trainingData = expand(os.path.join(outdir, "{sample}", "{sample}.AutoRTtrainingData.txt"), sample=sampleInfo_dict.keys()),
        DeepLC_trainingData = expand(os.path.join(outdir, "{sample}", "{sample}.DeepLCtrainingData.csv"), sample=sampleInfo_dict.keys()),
        filterTraining_stat = expand(os.path.join(outdir, "{sample}", "{sample}.filterTrainingDataNum.stat.txt"), sample=sampleInfo_dict.keys()),
        AutoRT_model = expand(os.path.join(outdir, "{sample}", "AutoRTtraingmodel", "model.json"), sample=sampleInfo_dict.keys()),
        AutoRT_RT_diff_result = expand(os.path.join(outdir, "{sample}", "{sample}.AutoRT.obsVSpredRTdiff.txt"), sample=sampleInfo_dict.keys()),
        AutoRT_concatRAW_result = expand(os.path.join(outdir, "{sample}", "{sample}.AutoRTobsVSpredRTdiff.RawPSM.txt"), sample=sampleInfo_dict.keys()),
        DeepLC_RT_diff_result = expand(os.path.join(outdir, "{sample}", "{sample}.DeepLC.obsVSpredRTdiff.txt"), sample=sampleInfo_dict.keys()),
        DeepLC_concatRAW_result = expand(os.path.join(outdir, "{sample}", "{sample}.DeepLCobsVSpredRTdiff.RawPSM.txt"), sample=sampleInfo_dict.keys()),
        sample_allRTsoft_result = expand(os.path.join(outdir, "{sample}", "{sample}.allRTsoft.RawPSM.txt"), sample=sampleInfo_dict.keys()),
        RT_aggregate_file = os.path.join(outdir, config["selfSet"]["prefix"] + ".addRetetionTime.PSM.txt"),
        ## all analysis result
        allAnalysis_result = os.path.join(outdir, config["selfSet"]["prefix"] + ".addRT.addSpecSim.PSM.txt")

#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
""" spectra similarity prediction snakemake pepline """
## create raw format mgf
rule rawMGFformat:
    input:
        psmFile = lambda wildcards: sampleInfo_dict[wildcards.sample]['psmFile']
    output:
        raw_format_mgf = os.path.join(outdir, "{sample}", "{sample}.raw.format.mgf")
    params:
        python = config["bindir"]["base_python"],
        rawMGFdir = config["selfSet"]["rawMGFdir"],
        py = config["bindir"]["formatRAWmgf"]
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.rawMGFformat.txt")

    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.py} --mgfRaw {params.rawMGFdir} --infile {input.psmFile} --outfile {output.raw_format_mgf}
        '''

## create pDeep input file
rule pDeepInputFile:
    input:
        psmFile = lambda wildcards: sampleInfo_dict[wildcards.sample]['psmFile']
    output:
        pDeep2_input = os.path.join(outdir, "{sample}", "{sample}.pDeep2_input.txt")
    params:
        python = config["bindir"]["base_python"],
        py = config["bindir"]["createpDeepInput_py"]
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.pDeepInputFile.txt")

    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.py} --infile {input.psmFile} --outfile {output.pDeep2_input}
        '''

## run pDeep2 prediction
rule pDeep2Prediction:
    input:
        pDeep2_input = os.path.join(outdir, "{sample}", "{sample}.pDeep2_input.txt")
    output:
        pDeep2_output = os.path.join(outdir, "{sample}", "{sample}.pDeep2_output.mgf")
    params:
        python = config["bindir"]["pDeep2_python"],
        py = config["bindir"]["pDeep2_predict_py"]
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.pDeep2Prediction.txt")

    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.py} -e 0.3 -i Lumos -in {input.pDeep2_input} -out {output.pDeep2_output}
        '''

## format pDeep2 result
rule pDeep2MGFformat:
    input:
        pDeep2_output = os.path.join(outdir, "{sample}", "{sample}.pDeep2_output.mgf")
    output:
        pDeep2_format_mgf = os.path.join(outdir, "{sample}", "{sample}.pDeep2_output.format.mgf")
    params:
        python = config["bindir"]["base_python"],
        py = config["bindir"]["formatpDeep2mgf"]
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.pDeep2MGFformat.txt")

    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.py} -i {input.pDeep2_output} -o {output.pDeep2_format_mgf}
        '''

## spectra similarity Score [raw vs pDeep]
rule SpectraSimilarity:
    input:
        raw_format_mgf = os.path.join(outdir, "{sample}", "{sample}.raw.format.mgf"),
        pDeep2_format_mgf = os.path.join(outdir, "{sample}", "{sample}.pDeep2_output.format.mgf"),
        psmFile = lambda wildcards: sampleInfo_dict[wildcards.sample]['psmFile']
    output:
        similarity_result = os.path.join(outdir, "{sample}", "{sample}.rawVSpDeep.spectra.similarity.txt"),
        sample_psm_result = os.path.join(outdir, "{sample}", "{sample}.addSpectraSimilarityScore.PSM.txt")
    params:
        python = config["bindir"]["similarity_python"],
        python_base = config["bindir"]["base_python"],
        py = config["bindir"]["ModifiedCosine_py"],
        addSpectraSimilarityScoretoPSM_py = config["bindir"]["addSpectraSimilarityScoretoPSM_py"]
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.SpectraSimilarity.txt")

    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.py} --raw {input.raw_format_mgf} --predict {input.pDeep2_format_mgf} --outfile {output.similarity_result}
        {params.python_base} {params.addSpectraSimilarityScoretoPSM_py} --infile {input.psmFile} --score {output.similarity_result} --outfile {output.sample_psm_result}
        '''   

rule SpecSim_aggregate:
    input:
        sample_psm_result = expand(os.path.join(outdir, "{sample}", "{sample}.addSpectraSimilarityScore.PSM.txt"), sample=sampleInfo_dict.keys())
    output:
        SpecSim_aggregate_file = os.path.join(outdir, config["selfSet"]["prefix"] + ".addSpectraSimilarityScore.PSM.txt")
    params:
        python = config["bindir"]["base_python"],
        concat_py = config['bindir']['concat_py']
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.concat_py} --infile '{input}' --outfile {output.SpecSim_aggregate_file}
        '''


#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
#####################################################################################################################################################################################################################
""" retention time prediction """
## create train and test Data work
rule RT_TrainTestData:
    input:
        psmFile = lambda wildcards: sampleInfo_dict[wildcards.sample]['psmFile']
    output:
        AutoRT_testData = os.path.join(outdir, "{sample}", "{sample}.AutoRTtestData.txt"),
        DeepLC_testData = os.path.join(outdir, "{sample}", "{sample}.DeepLCpredData.csv"),
        AutoRT_trainingData = os.path.join(outdir, "{sample}", "{sample}.AutoRTtrainingData.txt"),
        DeepLC_trainingData = os.path.join(outdir, "{sample}", "{sample}.DeepLCtrainingData.csv"),
        filterTraining_stat = os.path.join(outdir, "{sample}", "{sample}.filterTrainingDataNum.stat.txt")
    params:
        python = config["bindir"]["base_python"],
        py = config["bindir"]["rtData_py"],
        outpath = os.path.join(outdir, "{sample}")
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.RT_TrainTestData.txt")
    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.py}   -i {input.psmFile}   -p {wildcards.sample}   -o {params.outpath}
        '''

## AutoRT traing
rule AutoRTtraining:
    input:
        AutoRT_trainingData = os.path.join(outdir, "{sample}", "{sample}.AutoRTtrainingData.txt"),
    output:
        AutoRT_model = os.path.join(outdir, "{sample}", "AutoRTtraingmodel", "model.json"),
    params:
        mplconfigdir = os.path.join(outdir, "{sample}", "MPLCONFIGDIR"),
        AutoRTtraingmodel = os.path.join(outdir, "{sample}", "AutoRTtraingmodel"),
        python = config["bindir"]["AutoRT_python"],
        py = config["bindir"]["AutoRT_py"],        
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.AutoRTtraining.txt")
    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        mkdir -p {params.mplconfigdir}
        mkdir -p {params.AutoRTtraingmodel}
        export MPLCONFIGDIR="{params.mplconfigdir}"
        {params.python}   {params.py}  train   -i {input.AutoRT_trainingData}   -o {params.AutoRTtraingmodel}   -m /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/xianghaitao/software/Proteomics/AutoRT/AutoRT-master/models/base_model/model.json   -e 100   -b 64   -u s   -sm min_max   -rlr   -n 20
        rm -rf {params.mplconfigdir}
        '''

## AutoRT predict RT
rule AutoRTtest:
    input:
        psmFile = lambda wildcards: sampleInfo_dict[wildcards.sample]['psmFile'],
        AutoRT_testData = os.path.join(outdir, "{sample}", "{sample}.AutoRTtestData.txt"),
        AutoRT_model = os.path.join(outdir, "{sample}", "AutoRTtraingmodel", "model.json"),
    output:
        AutoRT_RT_diff_result = os.path.join(outdir, "{sample}", "{sample}.AutoRT.obsVSpredRTdiff.txt"),
        AutoRT_concatRAW_result = os.path.join(outdir, "{sample}", "{sample}.AutoRTobsVSpredRTdiff.RawPSM.txt"),
    params:
        AutoRT_python = config["bindir"]["AutoRT_python"],
        base_python = config["bindir"]["base_python"],
        AutoRT_py = config["bindir"]["AutoRT_py"],
        AutoRT_RTdiff_py = config["bindir"]["AutoRT_RTdiff_py"],
        AutoRT_concatRAW_py = config["bindir"]["AutoRT_concatRAW_py"],
        outpath = os.path.join(outdir, "{sample}")
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.AutoRTtest.txt")
    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.AutoRT_python}   {params.AutoRT_py}   predict   -t {input.AutoRT_testData}   -s {input.AutoRT_model}   -o {params.outpath}   -p {wildcards.sample}.AutoRT.RTpred
        {params.base_python}   {params.AutoRT_RTdiff_py}   --infile {params.outpath}/{wildcards.sample}.AutoRT.RTpred.tsv   --outfile   {output.AutoRT_RT_diff_result}
        {params.base_python}   {params.AutoRT_concatRAW_py}   --rawpsm {input.psmFile}   --rtfile {output.AutoRT_RT_diff_result}   --outfile {output.AutoRT_concatRAW_result} --prefix AutoRT
        '''        


## DeepLC predict RT
rule DeepLC:
    input:
        psmFile = lambda wildcards: sampleInfo_dict[wildcards.sample]['psmFile'],
        DeepLC_trainingData = os.path.join(outdir, "{sample}", "{sample}.DeepLCtrainingData.csv"),
        DeepLC_testData = os.path.join(outdir, "{sample}", "{sample}.DeepLCpredData.csv"),
    output:
        DeepLC_RT_diff_result = os.path.join(outdir, "{sample}", "{sample}.DeepLC.obsVSpredRTdiff.txt"),
        DeepLC_concatRAW_result = os.path.join(outdir, "{sample}", "{sample}.DeepLCobsVSpredRTdiff.RawPSM.txt"),
    params:
        DeepLC_python = config["bindir"]["DeepLC_python"],
        base_python = config["bindir"]["base_python"],
        DeepLC_RTdiff_py = config["bindir"]["DeepLC_RTdiff_py"],
        DeepLC_concatRAW_py = config["bindir"]["DeepLC_concatRAW_py"],
        outpath = os.path.join(outdir, "{sample}")
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.DeepLC.txt")
    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.DeepLC_python}   --file_pred {input.DeepLC_testData}   --file_cal {input.DeepLC_trainingData}   --file_pred_out {params.outpath}/{wildcards.sample}.DeepLC.RTpred.csv --plot_predictions
        {params.base_python}   {params.DeepLC_RTdiff_py}   --infile {params.outpath}/{wildcards.sample}.DeepLC.RTpred.csv --outfile {output.DeepLC_RT_diff_result}
        {params.base_python}   {params.DeepLC_concatRAW_py}   --rawpsm {input.psmFile}   --rtfile {output.DeepLC_RT_diff_result}   --outfile {output.DeepLC_concatRAW_result} --prefix DeepLC
        '''

## merge AutoRT and DeepLC obsVSpredRTdiff.RawPSM > allRTsoft PSM
rule RTsoftMerge:
    input:
        AutoRT_concatRAW_result = os.path.join(outdir, "{sample}", "{sample}.AutoRTobsVSpredRTdiff.RawPSM.txt"),
        DeepLC_concatRAW_result = os.path.join(outdir, "{sample}", "{sample}.DeepLCobsVSpredRTdiff.RawPSM.txt"),
    output:
        sample_allRTsoft_result = os.path.join(outdir, "{sample}", "{sample}.allRTsoft.RawPSM.txt"),
    params:
        base_python = config["bindir"]["base_python"],
        mergeAutoRTDeepLC_py = config["bindir"]["mergeAutoRTDeepLC_py"],
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}.rule.RTsoftMerge.txt")
    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.base_python}   {params.mergeAutoRTDeepLC_py}   --AutoRT {input.AutoRT_concatRAW_result}   --DeepLC {input.DeepLC_concatRAW_result}   --outfile {output.sample_allRTsoft_result}
        '''    

rule RT_aggregate:
    input:
        sample_psm_result = expand(os.path.join(outdir, "{sample}", "{sample}.allRTsoft.RawPSM.txt"), sample=sampleInfo_dict.keys())
    output:
        RT_aggregate_file = os.path.join(outdir, config["selfSet"]["prefix"] + ".addRetetionTime.PSM.txt")
    params:
        python = config["bindir"]["base_python"],
        concat_py = config['bindir']['concat_py']
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.concat_py} --infile '{input}' --outfile {output.RT_aggregate_file}
        '''


rule RT_SpecSim_aggregate:
    input:
        SpecSim_aggregate_file = os.path.join(outdir, config["selfSet"]["prefix"] + ".addSpectraSimilarityScore.PSM.txt"),
        RT_aggregate_file = os.path.join(outdir, config["selfSet"]["prefix"] + ".addRetetionTime.PSM.txt")
    output:
        allAnalysis_result = os.path.join(outdir, config["selfSet"]["prefix"] + ".addRT.addSpecSim.PSM.txt")
    params:
        python = config["bindir"]["base_python"],
        concat_py = config['bindir']['RT_SpecSim_merge_py']
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.concat_py} --infile '{input}' --outfile {output.allAnalysis_result}
        '''