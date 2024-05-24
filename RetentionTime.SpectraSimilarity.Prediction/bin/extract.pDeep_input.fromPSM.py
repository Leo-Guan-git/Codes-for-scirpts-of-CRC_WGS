from optparse import OptionParser
import pandas as pd
import os, re, glob

'''Function'''
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

'''main'''
def main():
    infile = options.infile
    outfile = options.outfile

    Info_data = pd.read_csv(infile, sep='\t')
    Info_data = Info_data[['Sequence', 'Modifications', 'Charge']].drop_duplicates().copy()

    Info_data = formatModification(Info_data)
    Info_data = Info_data.rename(columns={'Sequence':'peptide', 'Modifications':'modification', 'Charge':'charge'})
    Info_data.to_csv(outfile, sep='\t', index=False)
    
    

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()