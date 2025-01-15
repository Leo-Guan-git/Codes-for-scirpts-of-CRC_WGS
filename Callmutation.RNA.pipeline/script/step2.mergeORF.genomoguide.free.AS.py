# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import numpy as np

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

input_info = pd.read_csv(sys.argv[1], sep='\t', header=0, encoding='utf-8')
# CRC01.AS.Tumor.pep8.Trinity.orfipy.merge.detail.info
info_df = pd.DataFrame(input_info)
mutAA = info_df.iloc[:, 2]
mutAAlength = info_df.iloc[:, 3]
mutSiteNum = info_df.iloc[:, 4]
mutSites = info_df.iloc[:, 5]
strand = info_df.iloc[:, 6]
geneNum = info_df.iloc[:, 7]
geneID = info_df.iloc[:, 8]
geneName = info_df.iloc[:, 9]
mutTypeNum = info_df.iloc[:, 10]
mutTypes = info_df.iloc[:, 11]
ORFInfo = info_df.iloc[:, 13]

input_pepid = pd.read_csv(sys.argv[2], sep='\t', header=None, encoding='utf-8')
# CRC01.Tumor.all.Trinity.orfipy.filtered.ID.pep
pepid_df = pd.DataFrame(input_pepid)

output = open(sys.argv[3], 'w')

output.write('MUTPEPid' + '\t' + 'mutAA' + '\t' + 'mutAAlength' + '\t' + 'mutSiteNum' + '\t' +
             'mutSites' + '\t' + 'strand' + '\t' + 'geneNum' + '\t' + 'geneID' + '\t' + 'geneName' + '\t' +
             'mutTypeNum' + '\t' + 'mutTypes' + '\t' + 'ORFInfo' + '\n')

for i in range(len(info_df)):
    MUTPEPid = pepid_df[pepid_df.iloc[:, 0].isin([mutAA[i]])]
    MUTPEPid1 = MUTPEPid.iloc[:, 1].tolist()
    MUTPEPid2 = str(MUTPEPid1).split("'")[1]
    output.write(MUTPEPid2 + '\t' + str(mutAA[i]) + '\t' + str(int(mutAAlength[i])) + '\t' + str(int(mutSiteNum[i])) +
                 '\t' + str(mutSites[i]) + '\t' + str(strand[i]) + '\t' + str(int(geneNum[i])) + '\t' + str(geneID[i]) +
                 '\t' + str(geneName[i]) + '\t' + str(int(mutTypeNum[i])) + '\t' + str(mutTypes[i]) +
                 '\t' + str(ORFInfo[i]) + '\n')

output.close()

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
