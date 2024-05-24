from optparse import OptionParser
import pandas as pd
import os, re, glob, math

def main():
    DeepLCresult_data = pd.read_csv(options.infile, index_col=0)
    # DeepLCresult_data = pd.read_csv(file, index_col=0)
    obsVSpredRTdiff_data = pd.DataFrame(
        {
            'Sequence': DeepLCresult_data['seq'],
            'obsRTinSEC': DeepLCresult_data['tr'],
            'Modifications': DeepLCresult_data['modifications'],
            'predRTinSEC': DeepLCresult_data['predicted_tr'],
            'deltaRTinSEC': [math.fabs(DeepLCresult_data['tr'][i]-DeepLCresult_data['predicted_tr'][i]) for i in range(len(DeepLCresult_data))]
        }
    )
    obsVSpredRTdiff_data.to_csv(options.outfile, sep = '\t', index = False)

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()
