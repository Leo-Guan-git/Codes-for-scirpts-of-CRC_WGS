from optparse import OptionParser
import pandas as pd
import os, re
from pyteomics import mgf

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
    return(mgf_df)

def main():
    infile = options.infile
    annofile = options.anno
    mgfFile = options.mgf
    outFile = options.outfile
    global  min_length, max_length
    min_length, max_length = options.aalength.split('-')[0],options.aalength.split('-')[1]
    Type = options.Type
    Sample = options.Sample

    # data pre-processing: [mask, Sequence]
    data = pd.read_csv(infile, sep = '\t', header=None)
    data = data[[0,5]]
    data = data.rename(columns = {0:'mask', 5:'Sequence'})
    data = data[data['Sequence'].apply(lambda seq : True if (len(seq)>=int(min_length) and len(seq)<=int(max_length)) else False)]
    data = data.drop_duplicates()
    data = data.reset_index(drop = True)
    
    # anno_data pre-processing: [Spectrum, mask]
    anno_data = pd.read_csv(annofile, sep = '\t')
    anno_data['Spectrum'] = anno_data['MS File'].map(str) + '.' + anno_data['ScanNum'].map(str)
    anno_data = anno_data.rename(columns = {'Prediction': 'mask'})
    anno_data = anno_data[['Spectrum', 'mask']]

    # data add Spectrum Info
    data = pd.merge(data, anno_data, how='left', on='mask')
    data = data[['Spectrum', 'Sequence']]
    data['Modifications'] = [','.join([str(loc.end()) + '|Oxidation' for loc in re.finditer('m',seq)]) if list(re.finditer('m', seq)) != [] else 'Unmodified' for seq in data['Sequence'].tolist()]
    data['Sequence'] = data['Sequence'].str.upper()
    data['Score'] = ['-' for i in range(len(data))]

    # get mgf file
    mgf_data = extract_MSMSInfo(mgfFile)

    PSM_data = pd.merge(data, mgf_data, how='inner', on='Spectrum')

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
    parser.add_option("-i", "--infile",dest="infile", help="SMSNet MS190169_146_m-mod_fdr5_against_SwissProt_AND_CRC01.tumor.mergedWithB.MUTPRO.tsv")
    parser.add_option("-a", "--anno",dest="anno", help="SMSNet MS190169_146_m-mod_fdr5.tsv")
    parser.add_option('-m', "--mgf", dest="mgf", help="give me a mgf to me")
    parser.add_option('-l', "--aalength", dest="aalength", help="the range length of amino acid")
    parser.add_option('-t', "--Type", dest="Type", help="the type of sample")
    parser.add_option('-s', "--Sample", dest="Sample", help="sample name")
    parser.add_option("-o", "--outfile", dest="outfile", help="give me a out file  to me")
    (options, args) = parser.parse_args()
    main()