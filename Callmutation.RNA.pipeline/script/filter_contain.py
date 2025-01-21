# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import numpy as np

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

input1 = pd.read_csv(sys.argv[1], sep='\t', header=None)
input1_df = pd.DataFrame(input1)

input2 = pd.read_csv(sys.argv[2], sep='\t', header=None)
input2_df = pd.DataFrame(input2)

output = open(sys.argv[3], 'w')

mutpeptide = input1_df.iloc[:, 0]

input2_df1 = input2_df[0].tolist()
input2_df2 = str(input2_df1)

for i in range(len(input1_df)):
    if str(mutpeptide[i]) not in input2_df2:
        output.write(str(mutpeptide[i]) + '\n')
        continue

output.close()

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))



