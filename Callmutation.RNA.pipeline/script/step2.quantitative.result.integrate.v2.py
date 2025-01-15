# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import argparse
import math

def main():
    # Read the data from the files
    input_Tumor = pd.read_csv(args.input1, sep='\t', header=0, encoding='utf-8')
    input_Normal = pd.read_csv(args.input2, sep='\t', header=0, encoding='utf-8')

    # Merge the dataframes based on common columns
    merged_data = input_Tumor.merge(input_Normal, on=['geneID', 'geneSymbol', 'geneType', 'annoEnsGeneID'], suffixes=('_Tumor', '_Normal'))

    # Calculate logFC directly on the merged dataframe
    merged_data['logFC'] = merged_data.apply(lambda row: "{:.2f}".format(math.log2(row['FPKM_Tumor'] / row['FPKM_Normal']), 2), axis=1)
    merged_data['FPKM_Tumor'] = merged_data['FPKM_Tumor'].apply(lambda x: "{:.2f}".format(x))
    merged_data['FPKM_Normal'] = merged_data['FPKM_Normal'].apply(lambda x: "{:.2f}".format(x))

    output_columns = ['sampleID_Tumor', 'geneID', 'geneSymbol', 'geneType', 'annoEnsGeneID', 'FPKM_Tumor', 'FPKM_Normal', 'logFC']
    output = merged_data[output_columns]

    # Create a copy of the DataFrame before renaming the columns
    output = output.copy()

    # Rename the columns of the output DataFrame
    output.rename(columns={				
        'sampleID_Tumor': 'sampleID',
        'FPKM_Tumor': 'tumorFPKM',
        'FPKM_Normal': 'normalFPKM',
    }, inplace=True)

    output.to_csv("{0}/{1}.expFeature.tsv".format(args.outdir,args.prefix), sep='\t', index=False)

if __name__ == "__main__":
    print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
    
    parser = argparse.ArgumentParser(description='Quantitative result integrate')
    parser.add_argument('--input1', type=str, help='Tumor Input File')
    parser.add_argument('--input2', type=str, help='Normal Input File')
    parser.add_argument('--prefix', type=str, help='Sample ID or name')
    parser.add_argument('--outdir', type=str, default='.', help='Output path')
    args = parser.parse_args()
    main()

    print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
