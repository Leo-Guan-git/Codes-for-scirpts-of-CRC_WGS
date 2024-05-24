from optparse import OptionParser
import pandas as pd
import os, re
from pyteomics import mgf

def get_row(infile):
    skipped_list = []
    with open(infile, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.split('\t')[0] == 'MTD':
                skipped_list.append(i)
    return (skipped_list) 

def extract_MSMSInfo(mgfFile):
    reader = mgf.read(mgfFile)
    spectrum_list = []
    rtinseconds_list = []
    pepmass_list = []
    charge_list = []
    for spectrum in reader:
        spectrum_list.append(spectrum['params']['title'].split(' ')[0].rsplit('.',2)[0])
        rtinseconds_list.append(spectrum['params']['rtinseconds'])
        pepmass_list.append(spectrum['params']['pepmass'][0])
        charge_list.append(str(spectrum['params']['charge'][0])[0])
    mgf_df = pd.DataFrame(
        {
            'Spectrum': spectrum_list,
            'RTsecond': rtinseconds_list,
            'Charge': charge_list,
            'Mass': pepmass_list
        }
    )
    mgf_df['Charge'] = mgf_df['Charge'].astype('int')
    return(mgf_df)

def stat_loc_mod(string, match_mass_list, df_modifications_mass):
    loc_mod_list = []
    tmp = 0
    for mass in match_mass_list:
        pattern_mass = re.sub('(.*)\+(.*)',r'\1[+]\2',mass)
        if re.search(pattern_mass,string) != None:
            loc_mod_list.append('|'.join([str(int(re.search(pattern_mass,string).start()) + 1 + tmp),df_modifications_mass.index[df_modifications_mass['mass'] == mass ][0]]))
            for m in loc_mod_list:
                if re.search(r'Acetyl|Carbamyl|Ammonia', m) != None:
                    tmp = -1
                    loc_mod_list.remove(m)
            string = string.replace(mass, 'n', 1)
    modifications = ','.join(loc_mod_list) if loc_mod_list != [] else 'Unmodified'
    return(modifications)

def change_Casanovo_modification(data):
    dict_modifications_mass = {
        'Carbamidomethyl': 'C+57.021',
        'Oxidation[M]': 'M+15.995',
        'Deamidated[N]': 'N+0.984',
        'Deamidated[Q]': 'Q+0.984',
        'Acetyl[ProteinN-term]': '+42.011',
        'Carbamyl[ProteinN-term]': '+43.006',
        'Ammonia-loss[N]': '-17.027'
    }
    df_modifications_mass = pd.DataFrame.from_dict(dict_modifications_mass,orient='index',columns=['mass'])
    keywords = [re.sub('(.*)\+(.*)',r'\1[+]\2',value) for value in list(dict_modifications_mass.values())]
    pattern = re.compile('|'.join(list(map(str,keywords))))
    for i in range(len(data)):
        if pattern.findall(data['Modifications'][i]) != []:
            match_mass_list = pattern.findall(data['Modifications'][i])
            modifications = stat_loc_mod(data['Modifications'][i], match_mass_list, df_modifications_mass)
            data.loc[i, 'Modifications'] = modifications
        else:
            data.loc[i, 'Modifications'] = 'Unmodified'
    return(data)

def mod_format(df):
    for i in range(len(df)):
        mod = df['Modifications'][i]
        mod = re.sub(r'Oxidation\[M\]', r'Oxidation', mod)
        mod = re.sub(r'Deamidated\[.\]', r'Deamidated', mod)
        mod = re.sub(r'Acetyl\[ProteinN-term\]', r'Acetyl', mod)
        mod = re.sub(r'Ammonia-loss\[N\]', r'Ammonia-loss', mod)
        mod = re.sub(r'Carbamyl\[ProteinN-term\]', r'Carbamyl', mod)
        df.loc[i, 'Modifications'] = mod
    return(df)

def read_CasanovoResult(infile, mgf_df):
    df = pd.read_csv(infile, sep = '\t', skiprows=get_row(infile))

    # removo PSM_ID = nan rows
    df = df.drop(df[df['PSM_ID'].isna()].index)

    # df PSM_ID == mgf_df spectrum index
    df['Spectrum'] = mgf_df.loc[df['PSM_ID'].tolist()]['Spectrum'].tolist()
    df['Modifications'] = df['sequence']    
    df['Score'] = df["search_engine_score[1]"]
    df['Charge'] = df['charge'].astype('int')
    df['Sequence'] = df['sequence'].apply(lambda seq: ''.join(re.findall(r'[A-Za-z]', str(seq))))
    df = df[df['Sequence'].apply(lambda seq : True if (len(seq)>=int(min_length) and len(seq)<=int(max_length)) else False)]# Sequence 8-15
    df = df[['Spectrum',  'Charge', 'Sequence', 'Score', 'Modifications']]
    df = df.reset_index(drop = True)
    # Modification format
    df = change_Casanovo_modification(df)
    df = mod_format(df)
    return(df)

def main():
    infile = options.infile
    mgfFile = options.mgf
    outFile = options.outfile
    global  min_length, max_length
    min_length, max_length = options.aalength.split('-')[0],options.aalength.split('-')[1]
    Type = options.Type
    Sample = options.Sample

    # get mgf file
    mgf_data = extract_MSMSInfo(mgfFile)
    
    data = read_CasanovoResult(infile, mgf_data)

    PSM_data = pd.merge(data, mgf_data, how='inner', on=['Spectrum', 'Charge'])
    PSM_data = PSM_data[['Spectrum', 'Charge', 'Mass', 'RTsecond', 'Sequence', 'Score', 'Modifications']]
    PSM_data.insert(loc = 0, column = 'Type', value = [Type for i in range(len(PSM_data))])
    PSM_data.insert(loc = 1, column = 'Sample', value = [Sample for i in range(len(PSM_data))])
    PSM_data.to_csv(os.path.abspath(outFile), sep = '\t', index = False)


if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="filiter PepNet result",
                            #version="1.1",
                            epilog="""\
*************
                            """)
    parser.add_option("-i", "--infile",dest="infile", help="MS2100930018_Casanovo_result.mztab")
    parser.add_option('-m', "--mgf", dest="mgf", help="give me a mgf to me")
    parser.add_option('-l', "--aalength", dest="aalength", help="the range length of amino acid")
    parser.add_option('-t', "--Type", dest="Type", help="the type of sample")
    parser.add_option('-s', "--Sample", dest="Sample", help="sample name")
    parser.add_option("-o", "--outfile", dest="outfile", help="give me a out file  to me")
    (options, args) = parser.parse_args()
    main()