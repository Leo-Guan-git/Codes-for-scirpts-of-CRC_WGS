import os
import pandas as pd
import datetime
"""
Function
"""
 
def create_pepline_workDir(infile, outdir):
     
    sampleInfo = pd.read_csv(infile, sep = '\t', header=None, names=['Type', 'Sample', 'mgffile', 'dbfile'])
    sampleInfo['mgf'] = sampleInfo['mgffile']  
    sampleInfo['SpectrumID'] = sampleInfo['mgf'].apply(lambda x: os.path.split(x)[1].split('.')[0])  
    sampleInfo['Sample_file_raw'] = [os.path.join(outdir, "raw", sampleInfo['Type'][i], sampleInfo['Sample'][i], sampleInfo['SpectrumID'][i]) for i in range(len(sampleInfo))]  
    sampleInfo['Sample_file_result'] = [os.path.join(outdir, "result", sampleInfo['Type'][i], sampleInfo['Sample'][i], sampleInfo['SpectrumID'][i]) for i in range(len(sampleInfo))]  
    sampleInfo['SampleMGF'] = [os.path.join(outdir, "raw", sampleInfo['Type'][i], sampleInfo['Sample'][i], sampleInfo['SpectrumID'][i], sampleInfo['SpectrumID'][i] + '.mgf') for i in range(len(sampleInfo))]  
    sampleInfo['SampleMGF_dir'] = sampleInfo['SampleMGF'].apply(lambda x: os.path.split(x)[0].rsplit('/',1)[0])
    sampleInfo['inpath'] = [os.path.join(outdir, "raw", sampleInfo['Type'][i], sampleInfo['Sample'][i]) for i in range(len(sampleInfo))]  
    sampleInfo['expname'] = sampleInfo['SpectrumID']  
    sampleInfo['database_path'] = sampleInfo['dbfile'].apply(lambda x: os.path.split(x)[0])  
    sampleInfo['database_name'] = sampleInfo['dbfile'].apply(lambda x: os.path.split(x)[1].rsplit('.',2)[0])  
    sampleInfo['ID'] = sampleInfo['Type'] + '#' +  sampleInfo['Sample'] + '#' + sampleInfo['SpectrumID']
    sampleInfo['SMSNet_result'] = outdir + '/' + sampleInfo['ID'] + '/' + config['SMSNet']['result'] + '/' + sampleInfo['SpectrumID']+'_m-mod_fdr5_against_'+sampleInfo['database_name']+'.tsv'
    sampleInfo_dict = sampleInfo.set_index("ID").to_dict(orient="index")    
    for i in range(len(sampleInfo)):
        if not os.path.exists(sampleInfo['Sample_file_raw'][i]):
            os.makedirs(sampleInfo['Sample_file_raw'][i])
         
         
        if not os.path.exists(sampleInfo['SampleMGF'][i]):
            os.system(f"ln -s {sampleInfo['mgf'][i]} {sampleInfo['SampleMGF'][i]}")
    return sampleInfo_dict, sampleInfo


"""
main
"""
infile = config['selfSet']['infile']
outdir = os.path.abspath(config["selfSet"]["outdir"])
 
 
 
sampleInfo_dict, sampleInfo = create_pepline_workDir(infile, outdir)
date = datetime.datetime.now().strftime('%Y%m%d')


"""
snakemake pepline
"""
rule all: 
    input: ##最终需要生成文件
        PepNet_denovo_result_file = expand(os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.prediction.tsv"), sample=sampleInfo_dict.keys()),
        PepNet_PSM_txt = expand(os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.PSM.txt"), sample=sampleInfo_dict.keys()),
        PepNet_PSM_addProtein = expand(os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.PSM.addProtein.txt"), sample=sampleInfo_dict.keys()),
        
        SMSNet_denovo_result_file = expand([os.path.join(outdir, "{sample}", config['SMSNet']['result'], "{Spectrum_ID}"+'_m-mod_fdr5_against_'+"{database_name}"+'.tsv')], zip, sample=sampleInfo_dict.keys(), Spectrum_ID=sampleInfo['SpectrumID'].to_list(), database_name=sampleInfo['database_name'].to_list()),
        SMSNet_PSM_txt = expand(os.path.join(outdir, "{sample}", config['SMSNet']['result'], "{sample}.SMSNet.PSM.txt"), sample=sampleInfo_dict.keys()),
        SMSNet_PSM_addProtein = expand(os.path.join(outdir, "{sample}", config['SMSNet']['result'], "{sample}.SMSNet.PSM.addProtein.txt"), sample=sampleInfo_dict.keys()),

        Casanovo_denovo_result_file = expand(os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}_Casanovo_result.mztab"), sample=sampleInfo_dict.keys()),
        Casanovo_PSM_txt = expand(os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}.Casanovo.PSM.txt"), sample=sampleInfo_dict.keys()),
        Casanovo_PSM_addProtein =  expand(os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}.Casanovo.PSM.addProtein.txt"), sample=sampleInfo_dict.keys()),

        aggregate_file = os.path.join(outdir, "allBatches.denovo.PSM."+date+".txt")

rule PepNet:
    input:
        mgf = lambda wildcards: sampleInfo_dict[wildcards.sample]['mgf']
    output:
        PepNet_denovo_result_file = os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.prediction.tsv")
    params:
        python = config["PepNet"]["PepNet_python"],
        PepNet_denovo = config['PepNet']['PepNet_denovo'],
        PepNet_model = config['PepNet']['PepNet_model']
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}/{sample}.PepNet.txt")
    log:
        os.path.join(outdir, "{sample}", "log/{sample}/{sample}.PepNet.log")
    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.PepNet_denovo} --input {input.mgf} --model {params.PepNet_model} --output {output.PepNet_denovo_result_file}
        '''

rule PepNet_post:
    input:
        PepNet_denovo_result_file = os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.prediction.tsv"),
        mgf = lambda wildcards: sampleInfo_dict[wildcards.sample]['mgf']
    output:
        PepNet_PSM_txt = os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.PSM.txt"),
        PepNet_PSM_addProtein = os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.PSM.addProtein.txt")
    params:
        python = config["bindir"]["base_python"],
        postpy = config['PepNet']['PepNet_postpy'],
        animo_length = config['selfSet']['animo_length'],
        Type = lambda wildcards: sampleInfo_dict[wildcards.sample]['Type'],
        Sample = lambda wildcards: sampleInfo_dict[wildcards.sample]['Sample'],
        addProtein_py = config['bindir']['addProtein_py'],
        reference = config['selfSet']['reference']
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}/{sample}.PepNet_post.txt")
    log:
        os.path.join(outdir, "{sample}", "log/{sample}/{sample}.PepNet_post.log")
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.postpy} --infile {input.PepNet_denovo_result_file} --mgf {input.mgf} --aalength {params.animo_length} --Type {params.Type} --Sample {params.Sample} --outfile {output.PepNet_PSM_txt}
        {params.python} {params.addProtein_py} --infile {output.PepNet_PSM_txt}  --ref {params.reference} --outfile {output.PepNet_PSM_addProtein}
        '''

rule SMSNet:
    input:
        SampleMGF = lambda wildcards: sampleInfo_dict[wildcards.sample]['SampleMGF']
    output:
        SMSNet_denovo_result_file = os.path.join(outdir, "{sample}", config['SMSNet']['result'], "{Spectrum_ID}"+'_m-mod_fdr5_against_'+"{database_name}"+'.tsv')
    params:
        python = config["SMSNet"]["SMSNet_python"],
        SMSNet_denovo = config['SMSNet']['SMSNet_denovo'],
        SMSNet_model = config['SMSNet']['SMSNet_model'],
        SMSNet_dbsearch = config['SMSNet']['SMSNet_dbsearch'],
        SMSNet_thread = config['SMSNet']['SMSNet_thread'],
        inpath = lambda wildcards: sampleInfo_dict[wildcards.sample]['SampleMGF_dir'],
        expname = lambda wildcards: sampleInfo_dict[wildcards.sample]['expname'],
        database_path = lambda wildcards: sampleInfo_dict[wildcards.sample]['database_path'],
        database_name = lambda wildcards: sampleInfo_dict[wildcards.sample]['database_name'],
        outpath = os.path.join(outdir, "{sample}", config['SMSNet']['result'])
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}/{sample}.{Spectrum_ID}.{database_name}.SMSNet.txt")
    log:
        os.path.join(outdir, "{sample}", "log/{sample}/{sample}.{Spectrum_ID}.{database_name}.SMSNet.log")
    threads: 1
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        cd {params.outpath}
        rm -rf ./log_ablation
        {params.python} {params.SMSNet_denovo} --model_dir {params.SMSNet_model} --inference_input_file {input.SampleMGF} --rescore
        {params.python} {params.SMSNet_dbsearch} --inpath {params.inpath} --expname {params.expname} --dbpath {params.database_path} --dbname {params.database_name} --outpath {params.outpath} --thread {params.SMSNet_thread}
        '''

rule SMSNet_post:
    input:
        SampleMGF = lambda wildcards: sampleInfo_dict[wildcards.sample]['SampleMGF'],
        SMSNet_result = lambda wildcards: sampleInfo_dict[wildcards.sample]['SMSNet_result']
    output:
        SMSNet_PSM_txt = os.path.join(outdir, "{sample}", config['SMSNet']['result'], "{sample}.SMSNet.PSM.txt"),
        SMSNet_PSM_addProtein = os.path.join(outdir, "{sample}", config['SMSNet']['result'], "{sample}.SMSNet.PSM.addProtein.txt")
    params:
        python = config["bindir"]["base_python"],
        postpy = config['SMSNet']['SMSNet_postpy'],
        animo_length = config['selfSet']['animo_length'],
        inpath = lambda wildcards: sampleInfo_dict[wildcards.sample]['SampleMGF_dir'],
        expname = lambda wildcards: sampleInfo_dict[wildcards.sample]['expname'],
        Type = lambda wildcards: sampleInfo_dict[wildcards.sample]['Type'],
        Sample = lambda wildcards: sampleInfo_dict[wildcards.sample]['Sample'],
        addProtein_py = config['bindir']['addProtein_py'],
        reference = config['selfSet']['reference']
         
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}/{sample}.SMSNet_post.txt")
    log:
        os.path.join(outdir, "{sample}", "log/{sample}/{sample}.SMSNet_post.log")
    resources:
        tmpdir = config['selfSet']['temp']

    shell:
        '''
        {params.python} {params.postpy} --infile {input.SMSNet_result} --anno {params.inpath}/{params.expname}_m-mod_fdr5.tsv --mgf {input.SampleMGF} --aalength {params.animo_length} --Type {params.Type} --Sample {params.Sample} --outfile {output.SMSNet_PSM_txt}
        {params.python} {params.addProtein_py} --infile {output.SMSNet_PSM_txt}  --ref {params.reference} --outfile {output.SMSNet_PSM_addProtein}
        '''

rule Casanovo:
    input:
        SampleMGF = lambda wildcards: sampleInfo_dict[wildcards.sample]['SampleMGF']
    output:
        Casanovo_denovo_result_file = os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}_Casanovo_result.mztab")
    params:
        python = config["Casanovo"]["Casanovo_python"],
        Casanovo_config = config['Casanovo']['Casanovo_config'],
        Casanovo_model = config['Casanovo']['Casanovo_model'],   
    threads: 1
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}/{sample}.Casanovo.txt")
    log:
        os.path.join(outdir, "{sample}", "log/{sample}/{sample}.Casanovo.log")
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} --mode=denovo --model={params.Casanovo_model} --config={params.Casanovo_config} --peak_path={input.SampleMGF} --output={output.Casanovo_denovo_result_file}
        '''

rule Casanovo_post:
    input:
        Casanovo_denovo_result_file = os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}_Casanovo_result.mztab"),
        SampleMGF = lambda wildcards: sampleInfo_dict[wildcards.sample]['SampleMGF']
    output:
        Casanovo_PSM_txt = os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}.Casanovo.PSM.txt"),
        Casanovo_PSM_addProtein = os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}.Casanovo.PSM.addProtein.txt")
    params:
        python = config["bindir"]["base_python"],
        postpy = config['Casanovo']['Casanovo_postpy'],
        animo_length = config['selfSet']['animo_length'],
        Type = lambda wildcards: sampleInfo_dict[wildcards.sample]['Type'],
        Sample = lambda wildcards: sampleInfo_dict[wildcards.sample]['Sample'],
        addProtein_py = config['bindir']['addProtein_py'],
        reference = config['selfSet']['reference']
    benchmark:
        os.path.join(outdir, "{sample}", "benchmarks/{sample}/{sample}.Casanovo_post.txt")
    log:
        os.path.join(outdir, "{sample}", "log/{sample}/{sample}.Casanovo_post.log")
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.postpy} --infile {input.Casanovo_denovo_result_file} --mgf {input.SampleMGF} --aalength {params.animo_length} --Type {params.Type} --Sample {params.Sample} --outfile {output.Casanovo_PSM_txt}
        {params.python} {params.addProtein_py} --infile {output.Casanovo_PSM_txt}  --ref {params.reference} --outfile {output.Casanovo_PSM_addProtein}
        '''

rule aggregate:
    input:
        PepNet_PSM_addProtein = expand(os.path.join(outdir, "{sample}", config['PepNet']['result'], "{sample}.PepNet.PSM.addProtein.txt"), sample=sampleInfo_dict.keys()),
        SMSNet_PSM_addProtein = expand(os.path.join(outdir, "{sample}", config['SMSNet']['result'], "{sample}.SMSNet.PSM.addProtein.txt"), sample=sampleInfo_dict.keys()),
        Casanovo_PSM_addProtein =  expand(os.path.join(outdir, "{sample}", config['Casanovo']['result'], "{sample}.Casanovo.PSM.addProtein.txt"), sample=sampleInfo_dict.keys()),
    output:
        aggregate_file = os.path.join(outdir, "allBatches.denovo.PSM."+date+".txt")
    params:
        python = config["bindir"]["base_python"],
        concat_py = config['bindir']['concat_py']
    resources:
        tmpdir = config['selfSet']['temp']
    shell:
        '''
        {params.python} {params.concat_py} --infile '{input}' --outfile {output.aggregate_file}
        '''

