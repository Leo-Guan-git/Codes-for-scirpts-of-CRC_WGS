from optparse import OptionParser
import pandas as pd
import os, re, glob

'''Function'''
'''main'''
def main():
    merge_data = pd.DataFrame()
    for f in infile.split(' '):
        data = pd.read_csv(f, sep='\t')
        if merge_data.shape[0] == 0:
            merge_data = data
            continue
        merge_data = pd.merge(merge_data, data, on=['SampleID', 'SampleType', 'Omics', 'Software', 'Spectrum', 'Mass', 'Charge', 'RTsecond', 'Modifications', 'Sequence', 'Length', 'Score', 'Qvalue', 'DeltaScore', 'Proteins', 'pepClass'], how='left')
    merge_data.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    infile = options.infile
    outfile = options.outfile
    main()