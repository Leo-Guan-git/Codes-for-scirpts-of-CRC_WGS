from optparse import OptionParser
import pandas as pd
import math


def main():
    AutoRT_output = pd.read_csv(options.infile, sep = '\t')
    obsVSpredRTdiff_data = pd.DataFrame(
        {
            'Sequence': AutoRT_output['x'],
            'obsRTinSEC': AutoRT_output['y'],
            'predRTinSEC': [y_pred*60 for y_pred in AutoRT_output['y_pred']],
            'deltaRTinSEC': [math.fabs(AutoRT_output['y'][i]-AutoRT_output['y_pred'][i]*60) for i in range(len(AutoRT_output))]
        }
    )
    obsVSpredRTdiff_data.to_csv(options.outfile, sep = '\t', index = False)

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="AutoRT reuslt: RNA.Tumor.CRC01.dbsearch_MaxQuant.RTpred.tsv")
    parser.add_option("-o", "--outfile", dest="outfile", help="myself format reuslt: RNA.Tumor.CRC01.dbsearch_MaxQuant.obsVSpredRTdiff.txt")
    (options, args) = parser.parse_args()
    main()