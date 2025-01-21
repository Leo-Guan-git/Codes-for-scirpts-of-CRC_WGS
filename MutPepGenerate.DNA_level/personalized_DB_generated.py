#! /usr/bin/env python
# from __future__ import unicode_literals
# -*- coding: utf-8 -*-
# @Author: Xiangyu Guan
# @Date:   2025-01-21
# @Last Modified by:   Xiangyu Guan
# @Last Modified time: 2025-01-21

import pandas as pd
import gzip
import os
import argparse
import re


def ReadFile(file, mode='rt'):
    if file.endswith('.gz'):
        f = gzip.open(file, mode)
    else:
        f = open(file, mode)
    return(f)


def WriteFile(file, mode='wt'):
    if file.endswith('.gz'):
        f = gzip.open(file, mode)
    else:
        f = open(file, mode)
    return(f)


def NucleExtract(CHROM, genome_fai, genome_fasta):
    START_byte = genome_fai.loc[CHROM, 'OFFSET']
    Length_byte = genome_fai.loc[CHROM, 'LENGTH'] // genome_fai.loc[CHROM,
                                                                    'LINEBASES'] + genome_fai.loc[CHROM, 'LENGTH']
    genome_fasta.seek(START_byte, 0)
    nucle_seq = genome_fasta.read(Length_byte).replace("\n", "").upper()
    return(nucle_seq)


def VCFFileRead(filename, select_cols=[0, 1, 3, 4]):
    vcfs = pd.read_table(filename, header=None,
                         usecols=select_cols, comment='#')
    if len(select_cols) == 4:
        vcfs.columns = ['CHROM', 'POS', 'REF', 'ALT']
    # for lines that include multi-alternatives, just select first one
    vcfs['ALT'] = vcfs.ALT.apply(lambda x: x.split(',')[0])
    return(vcfs)


def GenomeFaiRead(filename):
    genome_fai = pd.read_table(filename, header=None)
    genome_fai.columns = ['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH']
    genome_fai.index = list(genome_fai['NAME'].values)
    return(genome_fai)


def NucleGenerate(raw_seq, variations):

    def Mutetype(x):
        if len(x.ALT) == len(x.REF):
            return('SNV')
        elif len(x.ALT) > len(x.REF):
            return('Insertion')
        else:
            return('Deletion')

    # ascending sort
    variations = variations.copy()
    variations = variations.sort_values(by='POS')
    # reindex
    variations.index = range(variations.shape[0])
    # type changing
    variations.POS = variations.POS.astype(int)
    # filter variations to make suer that all positions are included only once.
    if variations.shape[0] > 1:
        # start position confirm: start from 1
        start_pos = variations.iloc[0:-1,
                                    ].apply(lambda x: x.POS + len(x.REF), axis=1)
        start_pos.index = variations.index[1:]
        start_pos[variations.index[0]] = 1
        variations['START_POS'] = start_pos
        # variations.insert(variations.shape[1], 'START_POS', start_pos)
        while variations.query('START_POS > POS').shape[0]:
            variations.drop(variations.query(
                'START_POS > POS').index, inplace=True)
            if variations.shape[0] == 1:
                variations['START_POS'] = 1
                break
            start_pos = variations.iloc[0:-1,
                                        :].apply(lambda x: x.POS + len(x.REF), axis=1)
            start_pos.index = variations.index[1:]
            start_pos[variations.index[0]] = 1
            variations['START_POS'] = start_pos
        variations.index = range(variations.shape[0])
    else:
        variations['START_POS'] = 1
        # classify and select
        # select germline or somatic
    variations['ALT'] = variations.apply(
        lambda x: x.ALT_y if x.indicator_column == 'right_only' else x.ALT_x, axis=1)
    # 0 for germline, 1 for somatic
    variations['is_somatic'] = variations.indicator_column.apply(
        lambda x: 0 if x == 'right_only' else 1)
    # classify as snv, insertion or deletion
    variations['MUT_TYPE'] = variations.apply(Mutetype, axis=1)

    # new index for variation position in personal genome
    pos_change = variations.iloc[0:-
                                 1].apply(lambda x: len(x.ALT) - len(x.REF), axis=1)
    if pos_change.shape[0] > 1:
        pos_change = pos_change.cumsum()
        pos_change.index = variations.index[1:]
        pos_change[0] = 0
        variations['NEW_POS'] = pos_change + variations.POS
    else:
        variations['NEW_POS'] = variations.POS
    new_seq = ''.join(variations.apply(
        lambda x: raw_seq[x.START_POS - 1:x.POS - 1] + x.ALT, axis=1).values)
    new_seq += raw_seq[(variations.POS.values[-1] +
                        len(variations.REF.values[-1]) - 1):]
    # all nucleotides which is no A, T, C, G, U, N will be replaced as N
    new_seq = re.sub("[^ATCGUN]", "N", new_seq, flags=re.I)
    return(variations.get(['CHROM', 'POS', 'REF', 'NEW_POS', 'ALT', 'is_somatic', 'MUT_TYPE']), new_seq)


def SeqGenerate(chrom_sequence, pos, pos_seq, pep_length):
    # 待完善
    # 突变位点位于染色体开头或结尾，前/后序列长度不足，无法截取？？
    ns_length = pep_length * 3
    start_pos = pos - ns_length + 1
    end_pos = pos + len(pos_seq) + ns_length - 2
    return(start_pos, end_pos, chrom_sequence[start_pos - 1:end_pos])


def RevComp(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A', 'N': 'N'}
    # global complement
    return(''.join([complement[i] for i in seq[::-1]]))


def PepGenerate(seq):
    global AA_condon
    # like Stop condons, triplet codons included 'N' will be translated to B
    return(''.join([AA_condon[n] if n in AA_condon else 'B' for n in [seq[i:i + 3].upper() for i in range(0, len(seq), 3) if i + 3 <= len(seq)]]))


def SixBoxTrans(start_pos, end_pos, seq, pep_length):
    aa_results = list()
    starts = list()
    strands = list()
    # for forward strand
    strand = '+'
    for i in range(3):
        aa_seq = PepGenerate(seq[i:])
        tmp_start = start_pos + i
        for n in range(len(aa_seq) - pep_length + 1):
            mut_pep = aa_seq[n:n + pep_length]
            if 'B' in mut_pep:
                continue
            ns_start = tmp_start + n * 3
            aa_results.append(mut_pep)
            starts.append(ns_start)
            strands.append(strand)
    # for reverse strand
    seq_rev = RevComp(seq)
    strand = '-'
    for i in range(3):
        aa_seq = PepGenerate(seq_rev[i:])
        tmp_start = end_pos - i
        for n in range(len(aa_seq) - pep_length + 1):
            mut_pep = aa_seq[n:n + pep_length]
            if 'B' in mut_pep:
                continue
            ns_start = tmp_start - n * 3 - pep_length * 3 + 1
            aa_results.append(mut_pep)
            starts.append(ns_start)
            strands.append(strand)
    return(strands, starts, aa_results)


def MutPeptideGenerate(new_seq, selected_variations, chrom_id, pep_length):
    # ns means nucleotide sequence
    global AA_condon
    # target_ns = selected_variations.apply(lambda x: SeqGenerate(
    #     new_seq, x.POS, x.ALT, pep_length), axis=1)
    target_ns = selected_variations.apply(lambda x: SeqGenerate(
        new_seq, x.NEW_POS, x.ALT, pep_length), axis=1)
    target_pep = target_ns.apply(
        lambda x: SixBoxTrans(x[0], x[1], x[2], pep_length))
    return(target_pep)


def PersonalGenomeGenrtate(sample_id, chrom_id, outdir, ref_seq, new_seq, selected_variations):
    global chroms
    head_in = True
    if chrom_id in chroms:
        write_mode = 'wt'
        title = chrom_id
    else:
        write_mode = 'at'
        title = 'other'

    out_file = os.path.join(
        outdir, "{}.{}.mut.fa".format(sample_id, title))
    if title == 'other' and os.path.exists(out_file):
        head_in = False
    out_f = WriteFile(out_file, write_mode)

    # wirte to file
    print("\t>>>{chrom_id} Origin Seq length: {ref_len:d}\n\t>>>{chrom_id} Personalized Seq Length: {new_len:d}".format(
        chrom_id=chrom_id, ref_len=len(ref_seq), new_len=len(new_seq)))
    out_f.write(">{}\n{}\n".format(chrom_id, "\n".join(
        [new_seq[i:i + 100] for i in range(0, len(new_seq), 100)])))
    out_f.close()
    selected_variations.to_csv(os.path.join(outdir, "{}.{}.mut.tsv.gz".format(sample_id, title)),
                               sep="\t", index=False, header=head_in, mode=write_mode)
    if not title == 'other':
        os.system("bgzip -@ 4 -f {}".format(out_file))
        os.system("samtools faidx {}.gz".format(out_file))


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--referene', dest='reference', type=str, required=False,
                        help='genome fasta file used.')
    parser.add_argument('-n', '--index', dest='index', type=str, required=False,
                        help='samtools index file of genome fasta file.')
    parser.add_argument('-i', '--input', dest='input', type=str, required=False,
                        help='haplotype vcf file.')
    parser.add_argument('-i2', '--input2', dest='input2', type=str, required=False,
                        help='snv vcf file.')
    parser.add_argument('-i3', '--input3', dest='input3', type=str, required=False,
                        help='indel vcf file.')
    parser.add_argument('-l', '--length', dest='length', type=str, nargs='+', required=False,
                        help='mutation peptide lengths to generate[default is 8-15].')
    parser.add_argument('-c', '--chr', dest='chr', type=str, nargs='+', required=False,
                        help='which chromosome[s] to be generated.')
    parser.add_argument('-s', '--sample', dest='sample', type=str, required=True,
                        help='sample name used in generated file\'s name.')
    parser.add_argument('-o', '--ourdir', dest='outDir', type=str, required=True,
                        help='write result file to this path.')

    args = parser.parse_args()
    # fasta reference parse
    if args.reference:
        genome_fa_file = args.reference
    else:
        if args.input:
            in_f = ReadFile(args.input)
            for line in in_f:
                if line.startswith('##reference='):
                    genome_fa_file = line.rstrip().split('file:')[-1]
                    break
            in_f.close()
    genome_fasta = ReadFile(genome_fa_file)

    # fasta reference samtools index file parse
    if args.index:
        genome_fai_file = args.index
    else:
        genome_fai_file = genome_fa_file + '.fai'
    # fai index file read in
    genome_fai = GenomeFaiRead(genome_fai_file)
    # all vcf variation files read in
    if args.input:
        GERMLINE = VCFFileRead(args.input)
    else:
        GERMLINE = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT'])
    if args.input2:
        SNVs = VCFFileRead(args.input2)
    else:
        SNVs = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT'])
    if args.input3:
        INDELs = VCFFileRead(args.input3)
    else:
        INDELs = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT'])
    # SNVs = VCFFileRead(args.input2)
    # INDELs = VCFFileRead(args.input3)
    SOMATIC = pd.concat([SNVs, INDELs], ignore_index=True)
    # all variations merging
    All_Changed_pos = pd.merge(SOMATIC, GERMLINE, how='outer', on=[
                               'CHROM', 'POS', 'REF'], indicator='indicator_column')
    if args.chr:
        CHROMs = args.chr
    else:
        CHROMs = All_Changed_pos.CHROM.unique()

    try:
        os.remove(os.path.join(
            args.outDir, "{}.other.mut.tsv.gz".format(args.sample)))
    except FileNotFoundError:
        pass

    try:
        os.remove(os.path.join(
            args.outDir, "{}.other.mut.fa".format(args.sample)))
    except FileNotFoundError:
        pass

    if args.length:
        pep_len = [int(i) for i in args.length]
    else:
        pep_len = range(8, 16)

    # generate personal sequence for each chromosome
    for Chrom in CHROMs:
        # skip all HLA, decoy and EBV sequence included in reference
        if Chrom.startswith('HLA') or Chrom.endswith('decoy') or Chrom == 'chrEBV':
            continue
        # get target chromosome sequence from reference
        selected_nucle_seq = NucleExtract(
            Chrom, genome_fai, genome_fasta)
        # get chromosome variations
        selected_variations = All_Changed_pos.query(
            'CHROM == "{}"'.format(Chrom))
        if selected_variations.shape[0] == 0:
            print("\t No variation founded in {}".format(Chrom))
            continue
        # make sure that each position occured only once in a chromosome
        selected_variations = selected_variations.drop_duplicates(subset=[
                                                                  'POS'])
        selected_variations, new_seq = NucleGenerate(
            selected_nucle_seq, selected_variations)
        # generate mutate peptides for Specific length
        for each in pep_len:
            results = MutPeptideGenerate(
                new_seq, selected_variations, Chrom, each)
            selected_variations['PEP_{:d}'.format(each)] = results
            print("\t>>>{:d} long mutate peptides in {} number: {:d}".format(
                each, Chrom, results.apply(lambda x: len(x[2])).sum()))
        PersonalGenomeGenrtate(args.sample, Chrom, args.outDir,
                               selected_nucle_seq, new_seq, selected_variations)

    genome_fasta.close()

    out_file = os.path.join(
        args.outDir, "{}.other.mut.fa".format(args.sample))
    if os.path.exists(out_file):
        os.system("bgzip -@ 4 -f {}".format(out_file))
        os.system("samtools faidx {}.gz".format(out_file))


if __name__ == '__main__':
    # B(reak) for stop condon
    AA_condon = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R',
                 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'TCT': 'S', 'TCC': 'S',
                 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S', 'ATT': 'I', 'ATC': 'I',
                 'ATA': 'I', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
                 'CTG': 'L', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GTT': 'V',
                 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
                 'ACG': 'T', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'AAT': 'N',
                 'AAC': 'N', 'GAT': 'D', 'GAC': 'D', 'TGT': 'C', 'TGC': 'C', 'CAA': 'Q',
                 'CAG': 'Q', 'GAA': 'E', 'GAG': 'E', 'CAT': 'H', 'CAC': 'H', 'AAA': 'K',
                 'AAG': 'K', 'TTT': 'F', 'TTC': 'F', 'TAT': 'Y', 'TAC': 'Y', 'ATG': 'M',
                 'TGG': 'W', 'TAG': 'B', 'TGA': 'B', 'TAA': 'B'
                 }

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
              'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
    main()

