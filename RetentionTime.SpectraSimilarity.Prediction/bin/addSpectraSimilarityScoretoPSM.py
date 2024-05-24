from optparse import OptionParser
import pandas as pd
import os, re, glob

'''Function'''
def read_score_file(score_file):
    df = pd.read_csv(score_file, sep='\t')
    df['Sequence'] = df['raw_spectral'].apply(lambda x: x.split('#')[0])
    df['Modifications'] = df['raw_spectral'].apply(lambda x: x.split('#')[1])
    df['Charge'] = df['raw_spectral'].apply(lambda x: int(x.split('#')[2]))
    df['Spectrum'] = df['raw_spectral'].apply(lambda x: x.split('#')[3])
    df = df.drop('pDeep2_pred', axis=1)
    return df

def formatModification(df):
    df = df.copy()
    for i in df.index:
        sequence = df['Sequence'][i]
        modification = df['Modifications'][i]
        format_modification = []
        if modification == 'Unmodified':
            format_modification.append('')
        else:
            for mod in modification.split(','):
                loc = mod.split('|')[0]
                m = mod.split('|')[1]
                amino_acid = sequence[int(loc)-1]
                if m == 'Carbamidomethyl':
                    format_modification.append(','.join([loc, 'Carbamidomethyl'+'['+amino_acid+']']))
                elif m == 'Deamidated':
                    format_modification.append(','.join([loc, 'Deamidated'+'['+amino_acid+']']))
                elif m == 'Oxidation':
                    format_modification.append(','.join([loc, 'Oxidation'+'['+amino_acid+']']))
        df.loc[i,'Modifications'] = ';'.join(format_modification)
    return df

def read_infile(infile):
    df = pd.read_csv(infile, sep='\t')
    df = formatModification(df)
    return df

'''main'''
def main():
    infile = options.infile
    score_file = options.score
    outfile = options.outfile

    score_data = read_score_file(score_file)
    psm_data = read_infile(infile)
    result = pd.merge(psm_data, score_data, how='left', on=['Sequence', 'Modifications', 'Charge', 'Spectrum']).drop('raw_spectral', axis=1)

    result['Modifications'] = pd.read_csv(infile, sep='\t')['Modifications']

    result.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="PSM.txt")
    parser.add_option("-s", "--score",dest="score", help="spectra similarity score.txt")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()