# coding=utf-8
# 和 make.detail.info.py 不同的地方是加了 ORF 的起始密码子和终止密码子信息统计
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

input_pepinfo = pd.read_csv(sys.argv[2], sep='\t', header=0, encoding='utf-8')
# CRC01.as.named.sigs.8.pep.less.info
pepinfo_df = pd.DataFrame(input_pepinfo)
mutpepID = pepinfo_df.iloc[:, 0]
mutSites = pepinfo_df.iloc[:, 2]
mutID = pepinfo_df.iloc[:, 3]

input_geneinfo = pd.read_csv(sys.argv[3], sep='\t', header=0, encoding='utf-8')
# CRC01.as.named.sigs
geneinfo_df = pd.DataFrame(input_geneinfo)
geneid = geneinfo_df.iloc[:, 1]
gene = geneinfo_df.iloc[:, 2]
chr = geneinfo_df.iloc[:, 3]
strand = geneinfo_df.iloc[:, 4]
muttype = geneinfo_df.iloc[:, 5]
mutSite = geneinfo_df.iloc[:, 14]

input_pep_ORFinfo = pd.read_csv(sys.argv[4], sep='\t', header=None, encoding='utf-8')
# CRC01.Tumor.pep8.delpeptermination.genomeGuideReads.orfipy.AA.test.result
pep_ORFinfo_df = pd.DataFrame(input_pep_ORFinfo)

input_ORFinfo = pd.read_csv(sys.argv[5], sep='\t', header=0, encoding='utf-8')
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

output = open(sys.argv[6], 'w')
# CRC01.AS.Tumor.pep8.Trinity.ORF.pep.detail.info.txt

## step2.整合突变肽和相关 ORF 的详细信息 ##
# 构建header
output.write('MUTPEPid' + '\t' + 'peptideID' + '\t' + 'mutAA' + '\t' + 'mutAAlength' + '\t' + 'mutSiteNum' + '\t' +
             'mutSites' + '\t' + 'strand' + '\t' + 'geneNum' + '\t' + 'geneID' + '\t' + 'geneName' + '\t' +
             'mutTypeNum' + '\t' + 'mutTypes' + '\t' + 'genomeType' + '\t' + 'ORFInfo' + '\n')

for i in range(len(input_df)):
    geneid1 = []
    gene1 = []
    chr1 = []
    strand1 = []
    muttype1 = []
    mutSite1 = []

    output.write('MUTPEP.AS.' + str(i+1) + '\t')

    ## 统计突变肽的详细信息 ##
    pepinfo_df1 = pepinfo_df[pepinfo_df.iloc[:, 1].isin([mutpep[i]])]
    mutpepID1 = pepinfo_df1.iloc[:, 0]
    mutpepID2 = []
    mutID1 = pepinfo_df1.iloc[:, 3]
    l_list = pepinfo_df1.index.tolist()
    for l in l_list:
        mutpepID2.append(mutpepID1[l])
        mutID2 = int(mutID1[l])
        chr1.append(chr[mutID2-1])
        mutSite1.append(mutSite[mutID2-1])
        strand1.append(strand[mutID2 - 1])
        geneid1.append(geneid[mutID2 - 1])
        gene1.append(gene[mutID2 - 1])
        muttype1.append(muttype[mutID2 - 1])

    mutpepID2_unique = list(set(mutpepID2))
    mutpepID2_unique1 = ','.join(str(i) for i in mutpepID2_unique)
    output.write(str(mutpepID2_unique1) + '\t')

    output.write(str(mutpep[i]) + '\t' + str(len(mutpep[i])) + '\t')

    mutSiteNum = len(pepinfo_df1.iloc[:, 2].unique())
    output.write(str(mutSiteNum) + '\t')

    chr1_unique = list(set(chr1))
    chr1_unique1 = ','.join(str(i) for i in chr1_unique)
    mutSite1_unique = list(set(mutSite1))
    mutSite1_unique1 = ','.join(str(i) for i in mutSite1_unique)
    output.write(str(chr1_unique1) + '.' + str(mutSite1_unique1) + '\t')

    strand1_unique = list(set(strand1))
    strand1_unique1 = ','.join(str(i) for i in strand1_unique)
    output.write(str(strand1_unique1) + '\t')

    gene1_unique = list(set(gene1))
    geneNum = len(gene1_unique)
    gene1_unique1 = ','.join(str(i) for i in gene1_unique)
    geneid1_unique = list(set(geneid1))
    geneid1_unique1 = ','.join(str(i) for i in geneid1_unique)
    output.write(str(geneNum) + '\t' + str(geneid1_unique1) + '\t' + str(gene1_unique1) + '\t')

    muttype1_unique = list(set(muttype1))
    muttypeNum = len(muttype1_unique)
    muttype1_unique1 = ','.join(str(i) for i in muttype1_unique)
    output.write(str(muttypeNum) + '\t' + str(muttype1_unique1) + '\t')

    output.write('\t')

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
