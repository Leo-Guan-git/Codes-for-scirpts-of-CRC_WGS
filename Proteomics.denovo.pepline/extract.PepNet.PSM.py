from optparse import OptionParser
import pandas as pd
import os
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

def strandard_PSM(data, mgf_data):
    df = data.copy()
    # Pass uselessInfo
    df = df[['TITLE','DENOVO','Score','PPM Difference']]
    df = df[(df['PPM Difference']>=-20) & (df['PPM Difference']<=20)]
    df = df.rename(columns = {'TITLE': 'Spectrum', 'DENOVO': 'Sequence'})

    # Pass Sequence = NA rows
    df = df.dropna(axis=0, how='any',subset='Sequence')
    # Sequence 8-15
    df = df[df['Sequence'].apply(lambda seq : True if (len(seq)>=int(min_length) and len(seq)<=int(max_length)) else False)]
    df = df.reset_index(drop = True)

    df['Spectrum'] = df['Spectrum'].apply(lambda x: x.split(' ')[0].rsplit('.',2)[0]).copy()
    df['Modifications'] = ['Unmodified' for i in range(len(df))]

    df = pd.merge(df, mgf_data, how='inner', on='Spectrum')

    df = df[['Spectrum', 'Charge', 'Mass', 'RTsecond', 'Sequence', 'Score', 'Modifications']]
    
    return(df)


def main():
    infile = options.infile
    mgfFile = options.mgf
    outFile = options.outfile
    global  min_length, max_length
    min_length, max_length = options.aalength.split('-')[0],options.aalength.split('-')[1]
    Type = options.Type
    Sample = options.Sample

    data=pd.read_csv(infile, sep = '\t')
    
    # get mgf file
    mgf_data = extract_MSMSInfo(mgfFile)

    PSM_data = strandard_PSM(data, mgf_data)
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
    parser.add_option("-i", "--infile",dest="infile", help="PepNet reuslt prediction.tsv")
    parser.add_option('-m', "--mgf", dest="mgf", help="give me a mgf to me")
    parser.add_option('-l', "--aalength", dest="aalength", help="the range length of amino acid")
    parser.add_option('-t', "--Type", dest="Type", help="the type of sample")
    parser.add_option('-s', "--Sample", dest="Sample", help="sample name")
    parser.add_option("-o", "--outfile", dest="outfile", help="give me a out file  to me")
    (options, args) = parser.parse_args()
    main()