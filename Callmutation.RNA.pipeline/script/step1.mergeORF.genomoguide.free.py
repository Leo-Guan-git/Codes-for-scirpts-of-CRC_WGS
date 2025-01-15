# coding=utf-8
import os
import sys
import datetime
import pandas as pd
import numpy as np

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))

input_guide_info = pd.read_csv(sys.argv[1], sep='\t', header=0, encoding='utf-8')
# CRC01.AS.Tumor.pep8.Trinity.genomeGuideReads.orfipy.detail.info
guide_df = pd.DataFrame(input_guide_info)

input_free_info = pd.read_csv(sys.argv[2], sep='\t', header=0, encoding='utf-8')
# CRC01.AS.Tumor.pep8.Trinity.STARunAligned.orfipy.detail.info
free_df = pd.DataFrame(input_free_info)

guide_df["ORFInfo"] = guide_df["ORFInfo"].map(str) + free_df["ORFInfo"].map(str)

guide_df.to_csv(sys.argv[3], sep='\t', header=True, index=False)
# CRC01.AS.Tumor.pep8.Trinity.orfipy.merge.detail.info

print(datetime.datetime.now().strftime('%Y-%m-%d  %H:%M:%S'))
