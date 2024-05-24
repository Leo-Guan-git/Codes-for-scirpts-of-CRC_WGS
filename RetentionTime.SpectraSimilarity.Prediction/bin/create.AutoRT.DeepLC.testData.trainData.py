from optparse import OptionParser
import pandas as pd
import os, re, glob

def createAutoRTtrainingData(x, y):
    AutoRTtraing_data = pd.DataFrame(
        {
            'x': [x],
            'y': [y]
        }
    )
    return(AutoRTtraing_data)

def transform_AutoRTestData(x, y):
    AutoRTtest_data = pd.DataFrame(
        {
            'x': x,
            'y': y
        }
    )
    return(AutoRTtest_data)

def createDeepLCtrainingData(seq, modifications, rt):
    DeepLCtraing_data = pd.DataFrame(
        {
            'seq': [seq],
            'modifications': [modifications],
            'tr': [rt]
        }
    )
    return(DeepLCtraing_data)

def transform_DeepLCTestData(seq, modifications, rt):
    DeepLCTest_data = pd.DataFrame(
        {
            'seq': seq,
            'modifications': modifications,
            'tr': rt
        }
    )
    return(DeepLCTest_data)

def gettrainPassPSM(sample_df):
    # PEPScore < 0.01
    df = sample_df[sample_df['PEPscore'] <= 0.01].copy()
    df = df.reset_index(drop = True)
    # removo only exist DENOVOPRO PSM
    save_index_list = []
    for i in range(len(df)):
        preteion_list = df['Proteins'][i].split(',')
        for p in preteion_list:
            if re.search('DENOVO', p) == None:
                save_index_list.append(i)
                break
    df = df.loc[save_index_list]
    df = df.reset_index(drop = True)
    
    # statistics Spectrum Sequence numbers
    stat_df = df[['Spectrum', 'Sequence', 'Modifications']].value_counts().reset_index()
    stat_df.rename(columns={0:'repetition_times'},inplace = True)
    
    allAutoRTtraining_data = pd.DataFrame()
    allDeepLCtraining_data = pd.DataFrame()
    
    for i in stat_df[stat_df['repetition_times'] == 3].index:
        spectrum = stat_df['Spectrum'][i]
        sequence = stat_df['Sequence'][i]
        modification = stat_df['Modifications'][i]
        rtinseconds = df[(df['Spectrum'] == spectrum) & (df['Sequence'] == sequence) & (df['Modifications'] == modification)]['RTsecond'].values[0]
        
        # AutoRT trainData
        AutoRTtraing_data = createAutoRTtrainingData(sequence, rtinseconds)
        allAutoRTtraining_data = pd.concat([allAutoRTtraining_data, AutoRTtraing_data])

        # DeepLC trainData
        DeepLCtraing_data = createDeepLCtrainingData(sequence, modification.replace(',','|').replace('Unmodified',''), rtinseconds)
        allDeepLCtraining_data = pd.concat([allDeepLCtraining_data, DeepLCtraing_data])    
    
    # trainData Uniq
    allAutoRTtraining_data = allAutoRTtraining_data.drop_duplicates()
    allDeepLCtraining_data = allDeepLCtraining_data.drop_duplicates()
    
    return(allAutoRTtraining_data, allDeepLCtraining_data) 


def main():
    global outdir, prefix
    sample_df = pd.read_csv(infile, sep = '\t')
    
    # create Test data
    AutoRTtest_data = transform_AutoRTestData(sample_df['Sequence'], sample_df['RTsecond'])
    DeepLCTest_data = transform_DeepLCTestData(sample_df['Sequence'], sample_df['Modifications'].str.replace(',','|').replace('Unmodified',''), sample_df['RTsecond'])
    AutoRTtest_data.to_csv(os.path.join(outdir, '.'.join([prefix, 'AutoRTtestData.txt'])), sep = '\t', index = False)
    DeepLCTest_data.to_csv(os.path.join(outdir, '.'.join([prefix, 'DeepLCpredData.csv'])), sep = ',', index = False)

    # create Training data
    AutoRTtraining_data, DeepLCtraining_data = gettrainPassPSM(sample_df)
    AutoRTtraining_data.to_csv(os.path.join(outdir, '.'.join([prefix, 'AutoRTtrainingData.txt'])), sep = '\t', index = False)
    DeepLCtraining_data.to_csv(os.path.join(outdir, '.'.join([prefix, 'DeepLCtrainingData.csv'])), sep = ',', index = False)

    traingDataStat_data = pd.DataFrame(
        {
            'AutoRTtrainDataNum': [len(AutoRTtraining_data)],
            'DeepLCtrainDataNum': [len(DeepLCtraining_data)]
        }
    )
    traingDataStat_data.to_csv(os.path.join(outdir, '.'.join([prefix, 'filterTrainingDataNum.stat.txt'])), sep = '\t', index = False)


if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-p", "--prefix", dest="prefix", help="")
    parser.add_option("-o", "--outdir", dest="outdir", help="")
    (options, args) = parser.parse_args()
    infile = options.infile
    prefix = options.prefix
    outdir = options.outdir
    main()