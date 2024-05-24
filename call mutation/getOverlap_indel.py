#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import gzip
import os
import sys
import re
import datetime

if not 4 < len(sys.argv) < 8:
    print("ERROR: Please Use as:")
    print("\"python {} [mutect2.vcf[.gz]] [strelka2.indels.vcf[.gz]] [strelka.indel.vcf[.gz]] [svaba.vcf[.gz]] merged.indel.vcf uniq.indel.vcf\"".format(sys.argv[0]))
    print("Please input at least [two] indel.vcf[.gz] files.")
    print("[merged.indel.vcf] and [uniq.indel.vcf] are two oupput files.")
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
        indel_num = 0
        header = f.readline().rstrip("\n").split("\t")
        for i in range(0, len(header)):
            if 'tumor' in header[i].lower():
                tumor_ind = i
            if 'normal' in header[i].lower():
                normal_ind = i
        if 'svaba' in vcf_file_name:
            for line in f:
                line = line.rstrip("\n").split("\t")
                CHROM, POS, REF, ALTS, FILTER, FORMAT = line[0], line[1], line[3].upper(
                ), line[4].upper(), line[6], line[8].split(":")
                TUMOR, NORMAL = line[tumor_ind].split(
                    ':'), line[normal_ind].split(':')
                if len(REF) != len(ALTS):
                    alt = int(TUMOR[FORMAT.index('AD')])
                    ref = int(TUMOR[FORMAT.index('DP')])
                    key = "{}\t{}\t.\t{}\t{}".format(CHROM, POS, REF, ALTS)
                    alt_freq = float(alt) / ref
                    Result_Update(ref, alt, alt_freq, key)
                    indel_num += 1
        elif 'strelka' in vcf_file_name:
            for line in f:
                line = line.rstrip("\n").split("\t")
                CHROM, POS, REF, ALTS, FILTER, FORMAT = line[0], line[1], line[3].upper(
                ), line[4].upper(), line[6], line[8].split(":")
                TUMOR, NORMAL = line[tumor_ind].split(
                    ':'), line[normal_ind].split(':')
                ref = int(TUMOR[FORMAT.index('TAR')].split(',')[0])
                alt = int(TUMOR[FORMAT.index('TIR')].split(',')[0])
                if len(REF) != len(ALTS):
                    key = "{}\t{}\t.\t{}\t{}".format(CHROM, POS, REF, ALTS)
                    # alt_freq = float(alt) / (ref + alt)
                    try:
                        alt_freq = float(alt) / (ref + alt)
                    except:
                        alt_freq = 0.0
                    Result_Update(ref, alt, alt_freq, key)
                    indel_num += 1
        elif 'mutect' in vcf_file_name:
            for line in f:
                line = line.rstrip("\n").split("\t")
                CHROM, POS, REF, ALTS, FILTER, FORMAT = line[0], line[1], line[3].upper(
                ), line[4].upper(), line[6], line[8].split(":")
                TUMOR, NORMAL = line[tumor_ind].split(
                    ':'), line[normal_ind].split(':')
                ALTS = ALTS.split(',')
                AF = [float(i) for i in TUMOR[FORMAT.index('AF')].split(',')]
                AD = [int(i) for i in TUMOR[FORMAT.index('AD')].split(',')]
                # tumor_GT = [int(i)
                #             for i in TUMOR[FORMAT.index('GT')].split('/')]
                # normal_GT = [int(i)
                #              for i in NORMAL[FORMAT.index('GT')].split('/')]
                tumor_GT = [int(i) for i in re.split(
                    '[/\|]', TUMOR[FORMAT.index('GT')])]
                normal_GT = [int(i) for i in re.split(
                    '[/\|]', NORMAL[FORMAT.index('GT')])]
                alt_GT = list(set(tumor_GT) - set(normal_GT))
                ref = sum(AD)
                for i in alt_GT:
                    if i != 0 and len(ALTS[i - 1]) != len(REF):
                        key = "{}\t{}\t.\t{}\t{}".format(
                            CHROM, POS, REF, ALTS[i - 1])
                        alt = AD[i]
                        alt_freq = float(AF[i - 1])
                        Result_Update(ref, alt, alt_freq, key)
                        indel_num += 1
        print("{} >>>> Get {} indels from {} <<<".format(
            datetime.datetime.now(), indel_num, os.path.basename(vcf_file)))
    vcf_handle.close()

f1 = open(merged_output, 'w')
f2 = open(uniq_output, 'w')
merged_indel, uniq_indel = 0, 0
print("{} >>>>>>> Output to {} <<<<<<<".format(
    datetime.datetime.now(), merged_output))
print("{} >>>>>>> Output to {} <<<<<<<".format(
    datetime.datetime.now(), uniq_output))
for key in sorted(stat.keys()):
    if stat[key] == 1:
        f2.write("{}\t.\tPASS\t{}\n".format(key, freq[key]))
        uniq_indel += 1
    else:
        f1.write("{}\t.\tPASS\t{}\n".format(key, freq[key]))
        merged_indel += 1
f1.close()
f2.close()
print("{} >>>> Get {} merged indels <<<".format(
    datetime.datetime.now(), merged_indel))
print("{} >>>> Get {} uniq indels <<<".format(
    datetime.datetime.now(), uniq_indel))
print("{} >>>>>>> Job Done ! <<<<<<<".format(datetime.datetime.now()))
