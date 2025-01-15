# coding=utf-8

import os
import sys
import datetime
import pandas as pd

start_time = datetime.datetime.now()


data = pd.read_csv(sys.argv[1], header=None, encoding='GBK', dtype=str, sep='\t')

row_num = len(data)

size = int(sys.argv[2])

j = 1

for start in range(0, row_num, size):
    stop = start + size
    filename = sys.argv[4].format(sys.argv[3], j)
    d = data[start: stop]
    print("Saving file : " + filename + ", data size : " + str(len(d)))
    d.to_csv(filename, encoding='GBK', header=None, index=None, sep='\t')
    j = j + 1

end_time = datetime.datetime.now()

print(start_time)
print(end_time)
