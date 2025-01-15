# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import numpy as np

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

input = pd.read_csv(sys.argv[1], sep='\t', header=None, encoding='utf-8')
input_df = pd.DataFrame(input)

mutpep = input_df.iloc[:, 0]

output = open(sys.argv[2], 'w')

for i in range(len(input_df)):
    output.write('>MUTPEP_Fusion_' + str(i+1) + '\n' + str(mutpep[i]) + '\n')

output.close()

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
