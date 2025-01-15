# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import argparse

def main():
    input = pd.read_csv(args.input, sep='\t', header=0, encoding='utf-8')
    gene_id = input.iloc[:, 0]
    FPKM = input.iloc[:, 6]

    input_gtf = pd.read_csv(args.gtf, sep='\t', header=None, encoding='utf-8')
    gene_info_gtf = input_gtf.iloc[:, 8]

    # Create a dictionary to store gene information from the GTF file
    gene_info_dict = {}
    for l in range(len(input_gtf)):
        gene_id_gtf = gene_info_gtf[l].split(';')[0]
        gene_id1_gtf = gene_id_gtf.split('"')[1]

        gene_type_gtf = gene_info_gtf[l].split(';')[1]
        gene_type1_gtf = gene_type_gtf.split('"')[1]

        gene_Symbol_gtf = gene_info_gtf[l].split(';')[2]
        gene_Symbol1_gtf = gene_Symbol_gtf.split('"')[1]

        gene_info_dict[gene_id1_gtf] = (gene_Symbol1_gtf, gene_type1_gtf)

    output = open("{0}/{1}.{2}.quantitative.integrate.txt".format(args.outdir,args.prefix,args.type), 'w')
    output.write('sampleID' + '\t' + 'geneID' + '\t' + 'geneSymbol' + '\t' + 'geneType' + '\t' + 'annoEnsGeneID' + '\t' + 'FPKM' + '\n')

    for i in range(len(input)):
        if FPKM[i] > 0:
            gene_id1 = str(gene_id[i])
            if gene_id1 in gene_info_dict:
                gene_Symbol2_gtf, gene_type2_gtf = gene_info_dict[gene_id1]

                output.write(str(args.prefix) + '\t' + str(gene_id[i]) + '\t' + str(gene_Symbol2_gtf) + '\t' + 
                             str(gene_type2_gtf) + '\t' + str(gene_id[i]) + '\t' + "{:.2f}".format(FPKM[i]) + '\n')

    output.close()

if __name__ == "__main__":
    print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
    
    parser = argparse.ArgumentParser(description='Quantitative result integrate')
    parser.add_argument('--input', type=str, help='Input File')
    parser.add_argument('--gtf', type=str, help='GTF File')
    parser.add_argument('--prefix', type=str, help='Sample ID or name')
    parser.add_argument('--type', type=str, help='Sample type')
    parser.add_argument('--outdir', type=str, default='.', help='Output path')
    args = parser.parse_args()
    main()

    print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
