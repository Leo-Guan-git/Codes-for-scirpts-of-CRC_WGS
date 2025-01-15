# coding=utf-8
import os
import sys
import datetime
import pandas as pd

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))


input = pd.read_csv(sys.argv[1], sep='\t', header=0)
# CRC01.AS.Tumor.all.Trinity.orfipy.merger.detail.info
input_df = pd.DataFrame(input)

MUTPEPid = input_df.iloc[:, 0]
mutAA = input_df.iloc[:, 1]
mutAAlength = input_df.iloc[:, 2]
mutSiteNum = input_df.iloc[:, 3]
mutSites = input_df.iloc[:, 4]
strand = input_df.iloc[:, 5]
geneNum = input_df.iloc[:, 6]
geneID = input_df.iloc[:, 7]
geneName = input_df.iloc[:, 8]
mutTypeNum = input_df.iloc[:, 9]
mutTypes = input_df.iloc[:, 10]
ORFInfo = input_df.iloc[:, 11]


output = open(sys.argv[2], 'w')
# CRC01.AS.Tumor.all.Trinity.orfipy.merger.detail.v2.info

output.write('MUTPEPid' + '\t' + 'mutAA' + '\t' + 'mutAAlength' + '\t' + 'mutSiteNum' + '\t' +
             'mutSites' + '\t' + 'strand' + '\t' + 'geneNum' + '\t' + 'geneID' + '\t' + 'geneName' + '\t' +
             'mutTypeNum' + '\t' + 'mutTypes' + '\t' + 'startCodes' + '\t' + 'ORFInfo' + '\n')

for i in range(len(input_df)):
    startCode1 = []
    # stopCode1 = []
    orf1 =  str(ORFInfo[i]).split('|')
    len_orf1 = len(orf1)
    for l in range(len_orf1):
        orf2 = orf1[l]
        if len(orf2) > 1:
            orf3 = str(orf2).split('.')
            startCode1.append(orf3[7])
            # stopCode1.append(orf3[8])
    
    startCode2 = '|'.join(set(startCode1))
    # stopCode1 = '|'.join(set(stopCode1))
    
    output.write(str(MUTPEPid[i]) + '\t' + str(mutAA[i]) + '\t' + str(int(mutAAlength[i])) + '\t' + str(int(mutSiteNum[i])) +
                 '\t' + str(mutSites[i]) + '\t' + str(strand[i]) + '\t' + str(int(geneNum[i])) + '\t' + str(geneID[i]) +
                 '\t' + str(geneName[i]) + '\t' + str(int(mutTypeNum[i])) + '\t' + str(mutTypes[i]) +
                 '\t' + str(startCode2) + '\t' + str(ORFInfo[i]) + '\n')

output.close()

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
