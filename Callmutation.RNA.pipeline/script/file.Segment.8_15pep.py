# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import numpy as np

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

input_pepinfo = pd.read_csv(sys.argv[1], sep='\t', header=None, encoding='utf-8')
# CRC01.Tumor.all.Trinity.orfipy.filtered.pep
pepinfo_df = pd.DataFrame(input_pepinfo)
mutpep = pepinfo_df.iloc[:, 0]

output8 = open(sys.argv[2], 'w')
output9 = open(sys.argv[3], 'w')
output10 = open(sys.argv[4], 'w')
output11 = open(sys.argv[5], 'w')
output12 = open(sys.argv[6], 'w')
output13 = open(sys.argv[7], 'w')
output14 = open(sys.argv[8], 'w')
output15 = open(sys.argv[9], 'w')

for i in range(len(pepinfo_df)):
    if len(mutpep[i]) == 8:
        output8.write(str(mutpep[i]) + '\n')
    elif len(mutpep[i]) == 9:
        output9.write(str(mutpep[i]) + '\n')
    elif len(mutpep[i]) == 10:
        output10.write(str(mutpep[i]) + '\n')
    elif len(mutpep[i]) == 11:
        output11.write(str(mutpep[i]) + '\n')
    elif len(mutpep[i]) == 12:
        output12.write(str(mutpep[i]) + '\n')
    elif len(mutpep[i]) == 13:
        output13.write(str(mutpep[i]) + '\n')
    elif len(mutpep[i]) == 14:
        output14.write(str(mutpep[i]) + '\n')
    elif len(mutpep[i]) == 15:
        output15.write(str(mutpep[i]) + '\n')

output8.close()
output9.close()
output10.close()
output11.close()
output12.close()
output13.close()
output14.close()
output15.close()

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
