# coding=utf-8
import os
import sys
import datetime
import re
import pandas as pd

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

input1 = pd.read_csv(sys.argv[1], sep='\t', header=None)
# CRC01.Fusion.Tumor.all.Trinity.orfipy.filtered.pep
input1_df = pd.DataFrame(input1)
    
mutpep1 = input1_df.iloc[:, 0]


output = open(sys.argv[2], 'w')
# CRC01.Fusion.Tumor.all.Trinity.orfipy.filtered.ID.pep

l = 1                                                                
for i in range(len(input1_df)):
    output.write(str(mutpep1[i]) + '\t' + 'MUTPEP.Fusion.' + str(l) + '\n' )
    l = l + 1


output.close()

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
