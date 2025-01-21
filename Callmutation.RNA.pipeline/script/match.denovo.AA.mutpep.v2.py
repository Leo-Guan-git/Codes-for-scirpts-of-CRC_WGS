# coding=utf-8
import os
import sys
import datetime
import re
import pandas as pd

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

info = pd.read_csv(sys.argv[1], sep='\t', header=None, encoding='utf-8')
#info = pd.read_csv('CRC01.merged.coorded.tsv.peps.9.info.test', sep='\t', header=0, encoding='utf-8')

input = pd.read_csv(sys.argv[2], sep='\t', header=None, encoding='utf-8')
# input = pd.read_csv('CRC01.StringTie.superreads.Tumor/Normal.ORFfinder.line.test.out', sep='\t', header=None, encoding='utf-8')
input_df = pd.DataFrame(input)

output = open(sys.argv[3], 'w')

input_doubleline = input_df[input_df.index % 2 == 1]

mutProtein = info.iloc[:, 0]
input_df1 = input_doubleline[0].tolist()
input_df2 = str(input_df1)

for i in range(len(info)):
    if str(mutProtein[i]) in input_df2:
        for l in re.finditer(str(mutProtein[i]), input_df2):
            str1 = input_df2[:l.end()]
            index1 = int((str1.count("'") + 1) / 2)
            output.write(str(mutProtein[i]) + '\t' + input_df.iloc[index1*2-2, :][0] + '\n')

output.close()
print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
