from optparse import OptionParser
import pandas as pd
import os, re, glob

def main():
    AutoRTfile = options.AutoRT
    DeepLCfile = options.DeepLC
    outfile = options.outfile

    AutoRT_data = pd.read_csv(AutoRTfile, sep = '\t')
    DeepLC_data = pd.read_csv(DeepLCfile, sep = '\t')
    
    merge_data = pd.merge(AutoRT_data, DeepLC_data)
    
    merge_data.to_csv(outfile, sep = '\t', index = False)

    
if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-a", "--AutoRT",dest="AutoRT", help="")
    parser.add_option("-d", "--DeepLC",dest="DeepLC", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()