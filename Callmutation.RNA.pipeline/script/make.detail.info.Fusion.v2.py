# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import numpy as np


print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

## step1.读取整合所需的数据表文件 ##
input = pd.read_csv(sys.argv[1], sep='\t', header=None, encoding='utf-8')
# CRC01.Tumor.all.Trinity.orfipy.filtered.pep
input_df = pd.DataFrame(input)
mutpep = input_df.iloc[:, 0]

input_pepinfo = pd.read_csv(sys.argv[2], sep='\t', header=None, encoding='utf-8')
# unfixrm_CRC01.Tumor..8.pep.less.info
pepinfo_df = pd.DataFrame(input_pepinfo)
# mutpepID = pepinfo_df.iloc[:, 0]
# strand = pepinfo_df.iloc[:, 2]
# mutgene = pepinfo_df.iloc[:, 3]
# mut_gene_info = pepinfo_df.iloc[:, 4]

# input_geneinfo = pd.read_csv(sys.argv[3], sep='\t', header=0, encoding='utf-8')
# CRC01.as.named.sigs
# geneinfo_df = pd.DataFrame(input_geneinfo)
# geneid = geneinfo_df.iloc[:, 1]
# gene = geneinfo_df.iloc[:, 2]
# chr = geneinfo_df.iloc[:, 3]
# strand = geneinfo_df.iloc[:, 4]
# muttype = geneinfo_df.iloc[:, 5]
# mutSite = geneinfo_df.iloc[:, 14]

input_pep_ORFinfo = pd.read_csv(sys.argv[3], sep='\t', header=None, encoding='utf-8')
# CRC01.Fusion.Tumor.pep8.genomeGuideReads.orfipy.AA.result
pep_ORFinfo_df = pd.DataFrame(input_pep_ORFinfo)

input_ORFinfo = pd.read_csv(sys.argv[4], sep='\t', header=0, encoding='utf-8')
# CRC01.Tumor.genomeGuideReads.ORF.detail.less.info
ORFinfo_df = pd.DataFrame(input_ORFinfo)
# transID = ORFinfo_df.iloc[:, 0]
# transOrfID = ORFinfo_df.iloc[:, 1]
# strand_orf = ORFinfo_df.iloc[:, 2]
# frame = ORFinfo_df.iloc[:, 3]
# start = ORFinfo_df.iloc[:, 4]
# end = ORFinfo_df.iloc[:, 5]
# length = ORFinfo_df.iloc[:, 6]
# startCode = ORFinfo_df.iloc[:, 7]
# stopCode = ORFinfo_df.iloc[:, 8]

output = open(sys.argv[5], 'w')
# CRC01.Fusion.Tumor.pep8.Trinity.ORF.pep.detail.info.txt

## step2.整合突变肽和相关 ORF 的详细信息 ##
# 构建header
output.write('MUTPEPid' + '\t' + 'peptideID' + '\t' + 'mutAA' + '\t' + 'mutAAlength' + '\t' + 'mutSiteNum' + '\t' +
             'mutSites' + '\t' + 'strand' + '\t' + 'geneNum' + '\t' + 'geneID' + '\t' + 'geneName' + '\t' + 'mutgene' + '\t'
             'mutTypeNum' + '\t' + 'mutTypes' + '\t' + 'ORFInfo' + '\n')

for i in range(len(input_df)):
    mutpepID2 = []
    geneid2 = []
    gene2 = []
    gene2_Num = []
    chr2 = []
    strand2 = []
    # muttype2 = []
    mutSite2 = []
    mutgene2 = []

    output.write('MUTPEP.Fusion.' + str(i+1) + '\t')

    ## 统计突变肽的详细信息 ##
    pepinfo_df1 = pepinfo_df[pepinfo_df.iloc[:, 1].isin([mutpep[i]])]
    mutpepID1 = pepinfo_df1.iloc[:, 0]
    strand1 = pepinfo_df1.iloc[:, 2]
    mutgene1 = pepinfo_df1.iloc[:, 3]
    mut_gene_info1 = pepinfo_df1.iloc[:, 4]
    
    l_list = pepinfo_df1.index.tolist()
    for l in l_list:
        mutpepID2.append(mutpepID1[l])
        mut_gene_info2 = str(mut_gene_info1[l]).split('>')[1]
        chr_info = str(mut_gene_info2).split(';')[2]
        chr1 = str(chr_info).split(':')[0]
        chr2.append(chr1)
        mutSite1 = str(chr_info).split(':')[1]
        mutSite2.append(mutSite1)
        strand2.append(strand1[l])
        gene1_info = str(mut_gene_info2).split(';')[1]
        gene2_info = str(mut_gene_info2).split(';')[3]
        gene_id1 = str(gene1_info).split('^')[1]
        gene_id2 = str(gene2_info).split('^')[1]
        gene_id3 = gene_id1,gene_id2
        gene_id4 = ','.join(gene_id3)
        geneid2.append(gene_id4)
        gene_1 = str(gene1_info).split('^')[0]
        gene2_Num.append(gene_1)
        gene_2 = str(gene2_info).split('^')[0]
        gene2_Num.append(gene_2)
        gene_3 = gene_1,gene_2
        gene_4 = ','.join(gene_3)
        gene2.append(gene_4)
        mutgene2.append(mutgene1[l])

    mutpepID2_unique = list(set(mutpepID2))
    mutpepID2_unique1 = ','.join(str(i) for i in mutpepID2_unique)
    output.write(str(mutpepID2_unique1) + '\t')

    output.write(str(mutpep[i]) + '\t' + str(len(mutpep[i])) + '\t')

    mutSite2_unique = list(set(mutSite2))
    mutSiteNum = len(mutSite2_unique)
    output.write(str(mutSiteNum) + '\t')

    mutSite2_unique1 = ','.join(str(i) for i in mutSite2_unique)
    chr2_unique = list(set(chr2))
    chr2_unique1 = ','.join(str(i) for i in chr2_unique)
    output.write(str(chr2_unique1) + '.' + str(mutSite2_unique1) + '\t')

    strand2_unique = list(set(strand2))
    strand2_unique1 = ','.join(str(i) for i in strand2_unique)
    output.write(str(strand2_unique1) + '\t')

    gene2_Num_unique = list(set(gene2_Num))
    geneNum = len(gene2_Num_unique)
    gene2_unique = list(set(gene2))
    gene2_unique1 = ','.join(str(i) for i in gene2_unique)
    geneid2_unique = list(set(geneid2))
    geneid2_unique1 = ','.join(str(i) for i in geneid2_unique)
    output.write(str(geneNum) + '\t' + str(geneid2_unique1) + '\t' + str(gene2_unique1) + '\t')

    mutgene2_unique = list(set(mutgene2))
    mutgene2_unique1 = ','.join(str(i) for i in mutgene2_unique)
    output.write(str(mutgene2_unique1) + '\t')

    # muttype1_unique = list(set(muttype1))
    # muttypeNum = len(muttype1_unique)
    # muttype1_unique1 = ','.join(str(i) for i in muttype1_unique)
    # output.write(str(muttypeNum) + '\t' + str(muttype1_unique1) + '\t')

    output.write('-' + '\t')

    output.write('-' + '\t')

    ## 统计突变肽对应的 ORF 详细信息 ##
    pep_ORFinfo_df1 = pep_ORFinfo_df[pep_ORFinfo_df.iloc[:, 0].isin([mutpep[i]])]
    mutpep_orf = pep_ORFinfo_df1.iloc[:, 0]
    orf = pep_ORFinfo_df1.iloc[:, 1]

    y_list = pep_ORFinfo_df1.index.tolist()
    for y in y_list:
        orf1 = str(orf[y]).split('>')[1]
        orf2 = str(orf1).split('_ORF')[0]
        ORFinfo_df1 = ORFinfo_df[ORFinfo_df.iloc[:, 0].isin([orf2])]
        transID1 = ORFinfo_df.iloc[:, 0]
        transOrfID1 = ORFinfo_df.iloc[:, 1]
        strand_orf1 = ORFinfo_df.iloc[:, 2]
        frame1 = ORFinfo_df.iloc[:, 3]
        start1 = ORFinfo_df.iloc[:, 4]
        end1 = ORFinfo_df.iloc[:, 5]
        length1 = ORFinfo_df.iloc[:, 6]
        startCode1 = ORFinfo_df.iloc[:, 7]
        stopCode1 = ORFinfo_df.iloc[:, 8]
        x_list = ORFinfo_df1.index.tolist()
        for x in x_list:
            output.write(str(transID1[x]) + '.' + str(transOrfID1[x]) + '.' + str(strand_orf1[x]) + '.' +
                         str(frame1[x]) + '.' + str(start1[x]) + '.' + str(end1[x]) + '.' + str(length1[x]) + '.' +
                         str(startCode1[x]) + '.' + str(stopCode1[x]) + '|')

    output.write('\n')

output.close()


print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
