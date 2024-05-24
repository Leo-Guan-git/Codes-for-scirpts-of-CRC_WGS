from optparse import OptionParser
import pandas as pd
import os, re, glob

def main():
    infile = options.infile
    outfile = options.outfile

    i = 1
    for f in infile.split(' '):
        print('start merge file: {}.'.format(f))
        data = pd.read_csv(f, sep = '\t')
        if i == 1:
            data.to_csv(outfile, sep = '\t', index = False)
            i += 1
        else:
            data.to_csv(outfile, sep = '\t', index = False, header = False, mode='a')

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()