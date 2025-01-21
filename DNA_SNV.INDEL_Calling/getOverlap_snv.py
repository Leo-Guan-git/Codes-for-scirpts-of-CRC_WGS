#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import gzip
import os
import sys
import re
import datetime

if not 4 < len(sys.argv) < 9:
    print("ERROR: Please Use as:")
    print("\"python {} [MuSE.vcf.gz] [mutect2.vcf.gz] [mutect.vcf.gz] [strelka2.snvs.vcf.gz] [strelka.snv.vcf] merged.snv.vcf uniq.snv.vcf\"".format(sys.argv[0]))
    print("Please input at least [two] snv.vcf files.")
    print("[merged.snv.vcf] and [uniq.snv.vcf] are two oupput files.")
    exit(1)

inputfiles = sys.argv[1:len(sys.argv) - 2]
merged_output = sys.argv[-2]
uniq_output = sys.argv[-1]

stat = dict()
freq = dict()


def Result_Update(ref, alt, alt_freq, key):
    if key in freq.keys():
        freq[key] += ";{},{},{:.3f}".format(ref, alt, alt_freq)
    else:
        freq[key] = "{},{},{:.3f}".format(ref, alt, alt_freq)
    if key in stat.keys():
        stat[key] += 1
    else:
        stat[key] = 1


for vcf_file in inputfiles:
    vcf_file_name = vcf_file.lower()
    vcf_handle = gzip.open(vcf_file) if vcf_file_name.endswith(
        ".gz") else open(vcf_file)
    print("{} >>>>>>> Reading {} <<<<<<<".format(
        datetime.datetime.now(), vcf_file))
    with vcf_handle as f:
        snv_num = 0
        header = f.readline()
        header = header.rstrip("\n").split("\t")
        for i in range(0, len(header)):
            if 'tumor' in header[i].lower():
                tumor_ind = i
            if 'normal' in header[i].lower():
                normal_ind = i
        if 'muse' in vcf_file_name:
            for line in f:
                line = line.rstrip("\n").split("\t")
                CHROM, POS, REF, ALTS, FILTER, FORMAT = line[0], line[1], line[3], line[4], line[6], line[8].split(
                    ":")
                TUMOR, NORMAL = line[tumor_ind].split(
                    ':'), line[normal_ind].split(':')
                ALTS = ALTS.split(',')
                AD = [int(i) for i in TUMOR[FORMAT.index('AD')].split(',')]
                ref = sum(AD)
                # tumor_GT = [int(i)
                #             for i in TUMOR[FORMAT.index('GT')].split('/')]
                # normal_GT = [int(i)
                #              for i in NORMAL[FORMAT.index('GT')].split('/')]
                tumor_GT = [int(i)
                            for i in re.split('[/\|]', TUMOR[FORMAT.index('GT')])]
                normal_GT = [int(i)
                             for i in re.split('[/\|]', NORMAL[FORMAT.index('GT')])]
                alt_GT = list(set(tumor_GT) - set(normal_GT))
                for i in alt_GT:
                    if i != 0 and len(ALTS[i - 1]) == len(REF):
                        key = "{}\t{}\t.\t{}\t{}".format(
                            CHROM, POS, REF, ALTS[i - 1])
                        alt = AD[i]
                        alt_freq = float(alt) / ref
                        Result_Update(ref, alt, alt_freq, key)
                        snv_num += 1
        elif 'strelka' in vcf_file_name:
                # AD = [int(TUMOR[FORMAT.index(i)].split(',')[0]) for i in ['AU','CU','GU','TU']]
                # ref = sum(AD)
            for line in f:
                line = line.rstrip("\n").split("\t")
                CHROM, POS, REF, ALTS, FILTER, FORMAT = line[0], line[1], line[3], line[4], line[6], line[8].split(
                    ":")
                TUMOR, NORMAL = line[tumor_ind].split(
                    ':'), line[normal_ind].split(':')
                ALTS = ALTS.split(',')
                AD = [int(TUMOR[FORMAT.index(i)].split(',')[0])
                      for i in ['AU', 'CU', 'GU', 'TU']]
                ref = sum(AD)
                for i in ALTS:
                    if i != 0 and len(i) == 1:
                        key = "{}\t{}\t.\t{}\t{}".format(CHROM, POS, REF, i)
                        # ref = int(TUMOR[FORMAT.index(REF + 'U')].split(',')[0])
                        alt = int(TUMOR[FORMAT.index(i + 'U')].split(',')[0])
                        # alt_freq = float(alt) / ref
                        try:
                            alt_freq = float(alt) / ref
                        except:
                            alt_freq = 0.0
                        Result_Update(ref, alt, alt_freq, key)
                        snv_num += 1
        elif 'mutect' in vcf_file_name:
            for line in f:
                line = line.rstrip("\n").split("\t")
                CHROM, POS, REF, ALTS, FILTER, FORMAT = line[0], line[1], line[3], line[4], line[6], line[8].split(
                    ":")
                TUMOR, NORMAL = line[tumor_ind].split(
                    ':'), line[normal_ind].split(':')
                ALTS = ALTS.split(',')
                AD = [int(i) for i in TUMOR[FORMAT.index('AD')].split(',')]
                if 'AF' in FORMAT:
                    AF = [float(i)
                          for i in TUMOR[FORMAT.index('AF')].split(',')]
                if 'FA' in FORMAT:
                    AF = [float(i)
                          for i in TUMOR[FORMAT.index('FA')].split(',')]
                # tumor_GT = [int(i)
                #             for i in TUMOR[FORMAT.index('GT')].split('/')]
                # normal_GT = [int(i)
                #              for i in NORMAL[FORMAT.index('GT')].split('/')]
                tumor_GT = [int(i)
                            for i in re.split('[/\|]', TUMOR[FORMAT.index('GT')])]
                normal_GT = [int(i)
                             for i in re.split('[/\|]', NORMAL[FORMAT.index('GT')])]
                alt_GT = list(set(tumor_GT) - set(normal_GT))
                ref = sum(AD)
                for i in alt_GT:
                    if i != 0 and len(ALTS[i - 1]) == len(REF):
                        tmp1 = []
                        for n in range(0, len(REF)):
                            if REF[n] != ALTS[i - 1][n]:
                                tmp1.append([i, n, REF[n], ALTS[i - 1][n]])
                        for m in range(0, len(tmp1)):
                            key = "{}\t{}\t.\t{}\t{}".format(CHROM, int(
                                POS) + tmp1[m][1], tmp1[m][2], tmp1[m][3])
                            alt = AD[i]
                            alt_freq = float(AF[i - 1])
                            Result_Update(ref, alt, alt_freq, key)
                            snv_num += 1
        print("{} >>>> Get {} snvs from {} <<<".format(
            datetime.datetime.now(), snv_num, os.path.basename(vcf_file)))
    vcf_handle.close()

f1 = open(merged_output, 'w')
f2 = open(uniq_output, 'w')
merged_snv, uniq_snv = 0, 0
print("{} >>>>>>> Output to {} <<<<<<<".format(
    datetime.datetime.now(), merged_output))
print("{} >>>>>>> Output to {} <<<<<<<".format(
    datetime.datetime.now(), uniq_output))
for key in sorted(stat.keys()):
    if stat[key] == 1:
        f2.write("{}\t.\tPASS\t{}\n".format(key, freq[key]))
        uniq_snv += 1
    else:
        f1.write("{}\t.\tPASS\t{}\n".format(key, freq[key]))
        merged_snv += 1
f1.close()
f2.close()
print("{} >>>> Get {} merged snvs <<<".format(
    datetime.datetime.now(), merged_snv))
print("{} >>>> Get {} uniq snvs <<<".format(datetime.datetime.now(), uniq_snv))
print("{} >>>>>>> Job Done ! <<<<<<<".format(datetime.datetime.now()))
