#!python
from optparse import OptionParser
import pandas as pd
import os, re, glob

def meregRawPSMandobsVSpredRTdiff(rawpsm, rtfile, prefix):
    rawpsm_df = pd.read_csv(rawpsm, sep = '\t')
    ref_df = pd.read_csv(rtfile, sep = '\t')
    # ref_df = ref_df.fillna('')
    rawpsm_df[prefix+'predRTinSEC'] = ref_df['predRTinSEC']
    rawpsm_df[prefix+'deltaRTinSEC'] = ref_df['deltaRTinSEC']
    # rawpsm_df_copy = rawpsm_df.copy()
    # # 1.17780 -> 1.1
    # rawpsm_df_copy['RTsecond'] = rawpsm_df_copy['RTsecond'].apply(lambda x: float('.'.join([str(x).split('.')[0],str(x).split('.')[1][0:1]])))
    # ref_df['obsRTinSEC'] = ref_df['obsRTinSEC'].apply(lambda x: float('.'.join([str(x).split('.')[0],str(x).split('.')[1][0:1]])))
    # for i in range(len(rawpsm_df_copy)):
    #     Sequence = rawpsm_df_copy['Sequence'][i]
    #     Rtinseconds = rawpsm_df_copy['RTsecond'][i]
    #     modification = rawpsm_df_copy['Modifications'][i].replace(',','|').replace('Unmodified','')
    #     predRTinSEC = ref_df[(ref_df['Sequence'] == Sequence) & (ref_df['obsRTinSEC'] == Rtinseconds) & (ref_df['Modifications'] == modification)]['predRTinSEC'].values[0]
    #     deltaRTinSEC = ref_df[(ref_df['Sequence'] == Sequence) & (ref_df['obsRTinSEC'] == Rtinseconds) & (ref_df['Modifications'] == modification)]['deltaRTinSEC'].values[0]
    #     rawpsm_df.loc[i, prefix+'predRTinSEC'] = predRTinSEC
    #     rawpsm_df.loc[i, prefix+'deltaRTinSEC'] = deltaRTinSEC
    return(rawpsm_df)        

    
def main():
    rawpsm = options.rawpsm
    rtfile = options.rtfile
    prefix = str(options.prefix)
    outfile = options.outfile

    merge_data = meregRawPSMandobsVSpredRTdiff(rawpsm, rtfile, prefix)
    merge_data.to_csv(outfile, sep = '\t', index = False)

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-r", "--rawpsm",dest="rawpsm", help="DNAmutDB.Tumor.CRC01.Comet.RawPSM.txt")
    parser.add_option("-t", "--rtfile",dest="rtfile", help="DNAmutDB.Tumor.CRC01.Comet.AutoRT.obsVSpredRTdiff.txt")
    parser.add_option("-n", "--prefix",dest="prefix", help="AutoRT or DeepLC")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()