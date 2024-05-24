from optparse import OptionParser
import pandas as pd
import os, re, glob

'''Function'''
'''main'''
def main():
    infile = options.infile
    outfile = options.outfile

    result = open(outfile, 'w')
    with open(infile,'r') as f:
        for line in f:
            line = line.strip()
            if line.split('=')[0] == 'TITLE':
                print('{}'.format(line.replace('|', '#')), file = result)
            elif line.split('=')[0] == 'pepinfo':
                continue
            elif line[0:2] == 'b+':
                continue
            elif line[0:2] == 'y+':
                continue
            elif line[0].isdigit() == True:
                print('{} {}'.format(line.split()[0], line.split()[1]),file = result)
            else:
                print('{}'.format(line), file = result)

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()