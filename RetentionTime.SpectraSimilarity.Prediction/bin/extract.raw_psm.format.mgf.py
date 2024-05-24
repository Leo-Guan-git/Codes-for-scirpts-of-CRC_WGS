from optparse import OptionParser
import pandas as pd
import os, re, glob
from pyteomics import mgf
import datetime

def read_mgf(mgf_file):
    Spertrum_list = []
    Sample_mgf = mgf.read(mgf_file)
    for spectrum in Sample_mgf:
        Spertrum = spectrum.get('params')["title"].split(' ')[0].rsplit('.',2)[0]
        Spertrum_list.append(Spertrum)
    return(Spertrum_list, Sample_mgf)

def formatModification(df):
    df = df.copy()
    for i in df.index:
        sequence = df['Sequence'][i]
        modification = df['Modifications'][i]
        format_modification = []
        if modification == 'Unmodified':
            format_modification.append('')
        else:
            for mod in modification.split(','):
                loc = mod.split('|')[0]
                m = mod.split('|')[1]
                amino_acid = sequence[int(loc)-1]
                if m == 'Carbamidomethyl':
                    format_modification.append(','.join([loc, 'Carbamidomethyl'+'['+amino_acid+']']))
                elif m == 'Deamidated':
                    format_modification.append(','.join([loc, 'Deamidated'+'['+amino_acid+']']))
                elif m == 'Oxidation':
                    format_modification.append(','.join([loc, 'Oxidation'+'['+amino_acid+']']))
        df.loc[i,'Modifications'] = ';'.join(format_modification)
    return df

def main():
    infile = options.infile
    mgfRaw = options.mgfRaw
    outfile = options.outfile
    
    data = pd.read_csv(infile, sep = '\t')    
    data = data[['Spectrum', 'Modifications', 'Sequence']].copy()
    data = data.sort_values(by='Spectrum')
    data = data.drop_duplicates()
    data = data.reset_index(drop=True)
    data = formatModification(data)

    result = open(outfile, 'w')
    for mgf_file in glob.glob(os.path.join(mgfRaw,'*')):
        Sample_id = mgf_file.split('/')[-1].split('.')[0]   # MS1800283_207
        if Sample_id not in data['Spectrum'].apply(lambda x: x.split('.')[0]).drop_duplicates().to_list():
            continue
        Spertrum_list, Sample_mgf = read_mgf(mgf_file)
        print('{} Begin time: {}'.format(Sample_id, datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        for df_index in data[data["Spectrum"].str.split('.',expand=True)[0]==Sample_id].index: # spectrum_id = MS1800283_207.9983
            spectrum_id = data['Spectrum'][df_index]
            modification =  data['Modifications'][df_index]
            sequence = data['Sequence'][df_index]
            Sample_mgf_index = Spertrum_list.index(spectrum_id)
            #Sample_mgf[Sample_mgf_index].get('params')['title'] = '.'.join([data.loc[data["Spectrum"] == spectrum_id]["Sequence"].tolist()[0], spectrum_id]) # SVFQTTDTK.MS1800283_207.9983
            Peptides_mgf = Sample_mgf[Sample_mgf_index]  # change MGF TiTLE 
            new_title = '#'.join([sequence, modification, str(Peptides_mgf.get('params')['charge'])[0], spectrum_id]) # SVFQTTDTK#MS1800283_207.9983#2
            Peptides_mgf.get('params')['title'] =  new_title
            #result = open(os.path.join(os.path.abspath(outdir), new_title+'.mgf'), 'w')
            # if len(list(Peptides_mgf.get('params')['pepmass'])) == 1:

            print(f'''BEGIN IONS
TITLE={new_title}
RTINSECONDS={Peptides_mgf.get('params')['rtinseconds']}
PEPMASS={' '.join([str(list(Peptides_mgf.get('params')['pepmass'])[0]), str(list(Peptides_mgf.get('params')['pepmass'])[1])]) if list(Peptides_mgf.get('params')['pepmass'])[1] != None else ' '.join([str(list(Peptides_mgf.get('params')['pepmass'])[0]), str(list(Peptides_mgf.get('params')['pepmass'])[0])])}
CHARGE={Peptides_mgf.get('params')['charge']}''', file = result)
            for i in range(len(Peptides_mgf.get('m/z array'))):
                print('{} {}'.format(Peptides_mgf.get('m/z array')[i], Peptides_mgf.get('intensity array')[i]), file = result)
            print('END IONS', file = result)
        print('{} Finished time: {}'.format(Sample_id, datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-m", "--mgfRaw",dest="mgfRaw", help="")
    parser.add_option("-i", "--infile",dest="infile", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()