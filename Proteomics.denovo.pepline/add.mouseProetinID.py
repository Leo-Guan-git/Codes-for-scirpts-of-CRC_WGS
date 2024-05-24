from optparse import OptionParser
import pandas as pd
import os, re
import gzip
from Bio import SeqIO

def add_proteins(fasta_file,data):
    handle = gzip.open(fasta_file, 'rt')
    dict_dataIndex_proteinID = {}
    for record in SeqIO.parse(handle, 'fasta'):
        fasta_proteinID = ':'.join([record.id.split('|')[1], record.id.split('|')[2].split('_')[0]])
        fasta_sequence = str(record.seq)
        for i in range(len(data)):
            self_sequence = data['Sequence'][i]
            if re.search(self_sequence.replace('I','L'),fasta_sequence.replace('I','L')) != None:
                dict_dataIndex_proteinID.setdefault(i,[]).append(fasta_proteinID)
    for i in range(len(data)):
        if i not in dict_dataIndex_proteinID:
            dict_dataIndex_proteinID[i] = '-'
        else:
            dict_dataIndex_proteinID[i] = ';'.join(dict_dataIndex_proteinID[i])
    data.insert(loc=len(data.columns), column='Proteins', value=[dict_dataIndex_proteinID[i] for i in range(len(data))])
    return(data)

def main():
    data = pd.read_csv(options.infile, sep = '\t', low_memory = False)
    data = add_proteins(options.ref, data)
    data.to_csv(options.outfile, sep = '\t', index = False)

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-r", "--ref",dest="ref", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()
