import sys,os
import tempfile, subprocess
from subprocess import PIPE, Popen
import pybedtools

def pyBedGetSeq(chr, start, end, strandness, myref):
    name = "name"
    score = "1"
    mybedline = chr + "\t" + start + "\t" + end + "\t" + "name" + "\t" + score + "\t" + strandness
    got_seq = "unk"
    myrefFa = myref
    myquery = pybedtools.BedTool(mybedline, from_string=True)
    got_seq = myquery.sequence(fi=myrefFa)
    try:
        myres = open(got_seq.seqfn).readlines()[1].strip()
    except Exception:
        myres = "nah"
    return  myres     #open(got_seq.seqfn).readlines()[1].strip()
# import pysam
table = """
TTT F      CTT L      ATT I      GTT V
TTC F      CTC L      ATC I      GTC V
TTA L      CTA L      ATA I      GTA V
TTG L      CTG L      ATG M      GTG V
TCT S      CCT P      ACT T      GCT A
TCC S      CCC P      ACC T      GCC A
TCA S      CCA P      ACA T      GCA A
TCG S      CCG P      ACG T      GCG A
TAT Y      CAT H      AAT N      GAT D
TAC Y      CAC H      AAC N      GAC D
TAA *      CAA Q      AAA K      GAA E
TAG *      CAG Q      AAG K      GAG E
TGT C      CGT R      AGT S      GGT G
TGC C      CGC R      AGC S      GGC G
TGA *      CGA R      AGA R      GGA G
TGG W      CGG R      AGG R      GGG G
"""

cwd = os.getcwd()

tableDict = dict(zip(table.split()[::2],table.split()[1::2]))

def getseq(chr, start, end, strandness, myref):
    name = "name"
    score = "1"
    mybedline = chr + "\t" + start + "\t" + end + "\t" + "name" + "\t" + score + "\t" + strandness
    got_seq = "unk"
    with  tempfile.NamedTemporaryFile(mode='w', dir= cwd) as fp:
        fp.write(mybedline)
        fp.seek(0)
        inputBed = fp.name
        syscmdline = "bedtools getfasta -s -tab -bed " + inputBed + " -fi " + myref
        cmdTask = subprocess.Popen(syscmdline,shell=True,stdout=subprocess.PIPE)
        mytext = cmdTask.stdout.readlines()
        got_seq = mytext[0].decode('UTF-8').split("\t")[1].strip()#.strip()
    return got_seq


def translate(dna):
    aa = ''
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if codon in tableDict.keys():
            aa += tableDict[codon]
        else:
            aa += "_"
    return aa

def reverse_complement(dna):
    revc = ""
    basepair = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    for c in dna:
        revc = basepair[c] + revc
    return revc

def six_frame_translate(dna, peptideLen):
    slideTimes = int(peptideLen)
    readStep = int(peptideLen) * 3
    revc = reverse_complement(dna)
    forward_amino_acids = []
    reverse_amino_acids = []
    for tranStart in range(readStep):
        forward_amino_acids.append(translate(dna[tranStart: tranStart+readStep]))
        reverse_amino_acids.append(translate(revc[tranStart: tranStart+readStep]))
    return (forward_amino_acids + reverse_amino_acids)

# FusionName	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint
# FBXO25--SEPTIN14	FBXO25^ENSG00000147364.17	chr8:435707:+	SEPTIN14^ENSG00000154997.9	chr7:55796092:-
# NRIP1--LINC02246	NRIP1^ENSG00000180530.11	chr21:15064745:-	LINC02246^ENSG00000281903.2	chr21:14857708:-

def ordering_translate_other(dna, peptideLen, strandness, leftLen):
    seqLen = len(dna)
    readStep = int(peptideLen) * 3
    lastSub = seqLen - int(peptideLen) * 3
    if strandness == "+":
        rcStrandness = "-"
    elif strandness == "-":
        rcStrandness = "+"
    else:
        strandness = "+"
        rcStrandness = "-"
    amino_acids = []
    relativeLoc = -(int(peptideLen) * 3 - 1)
    for subSeqStart in range(lastSub):
        finalBase = subSeqStart + readStep
        subSeq = dna[subSeqStart:finalBase]
        subPep = translate(subSeq)
        if subSeqStart >= int(leftLen):
            locTag = "on_right_gene"
        else:
            locTag = "on_left_gene"
        subInfo = subSeq + ":" + subPep + ":" + strandness + ":" + str(subSeqStart) + "-" + str(finalBase) + ":" + "reLoc=" + str(relativeLoc) + ":" + locTag
        amino_acids.append(subInfo)
        relativeLoc += 1
    # revDNA = reverse_complement(dna)
    rightLen = len(dna) - int(leftLen)
    relativeLocRev = int(peptideLen) * 3 - 1
    for subSeqStart in range(lastSub):
        finalBase = subSeqStart + readStep
        subSeq = dna[subSeqStart:finalBase]
        subPep = translate(subSeq)
        if subSeqStart < rightLen:
            locTag = "on_right_gene"
        else:
            locTag = "on_left_gene"
        subInfo = subSeq + ":" + subPep + ":" + rcStrandness + ":" + str(subSeqStart) + "-" + str(finalBase) + ":" + "reLoc=" + str(relativeLocRev) + ":" + locTag
        amino_acids.append(subInfo)
        relativeLocRev -= 1
    return amino_acids

def expandFusion(tumorSeq, strandness, pepLen, leftLen):
    outseqs = ordering_translate_other(tumorSeq, str(pepLen), strandness, leftLen )
    return outseqs


if __name__ == '__main__':
    myGeneDict = {}
    myGtfText = open(sys.argv[3]).readlines()
    for line in myGtfText:
        if "#" not in line:
            gtmp = line.split("\t")
            if gtmp[2] == "gene":
                gloci = gtmp[0] + ":" + gtmp[3] + ":" + gtmp[4] + ":" + gtmp[5]
                gID = gtmp[8].split()[1].lstrip("\"").strip(";").strip("\"")
                myGeneDict[gID] = gloci
    preMergedFile = open(sys.argv[1], 'r').readlines()[1:]
    # star-fusion.fusion_predictions.abridged.coding_effect.tsv
    pepLen = int(sys.argv[2])
    myrefFa = sys.argv[4]
    peptideID = 1
    eventID = 1
    readLen = pepLen * 3
    for line in  preMergedFile:
        tmp = line.split("\t")
        geneInfo =  tmp[0] + ";" + ";".join(tmp[6:9])
        leftGeneID = tmp[6].split("^")[1]
        rightGeneID = tmp[8].split("^")[1]
        leftGeneEND = tmp[7].split(":")[1]     ## Left  break location 
        rightGeneStart = tmp[9].split(":")[1]  ## Right  break location

        leftChr = tmp[7].split(":")[0]
        leftStart = myGeneDict[leftGeneID].split(":")[1]  ## GTF file, gene start
        leftStrand = tmp[7].split(":")[2]

        rightChr = tmp[9].split(":")[0]
        rightEND = myGeneDict[rightGeneID].split(":")[2]  ## GTF file, gene END  
        rightStrand  =  tmp[9].split(":")[2]

        if leftStrand == "+":
            leftSeq = pyBedGetSeq(leftChr, leftStart, leftGeneEND, leftStrand, myrefFa)
        else:
            leftStart = leftGeneEND
            leftGeneEND = myGeneDict[leftGeneID].split(":")[2]
            leftSeq = pyBedGetSeq(leftChr, leftStart, leftGeneEND, leftStrand, myrefFa)
        if rightStrand == "+":
            rightSeq = pyBedGetSeq(rightChr, rightGeneStart, rightEND, rightStrand, myrefFa)
        else:
            rightGeneStartMinus = myGeneDict[rightGeneID].split(":")[1]
            rightEND = rightGeneStart
            rightSeq = pyBedGetSeq(rightChr, rightGeneStartMinus, rightEND,  rightStrand, myrefFa)

        leftLen = len(leftSeq)
        rightLen = len(rightSeq)

        myTumorDNAsequence = leftSeq + rightSeq
        
        #deletedSeqLen = len(normalSeq) - 1
        #print("\t".join(["geneInfo", "leftGeneLength", "RightGeneLength", "FusionGeneSequence"]))
        #print(">" + geneInfo + "\t" + leftSeq + "\t" + rightSeq + "\t" + str(leftLen) + "\t" + str(rightLen), end = "\t")
        #print(myTumorDNAsequence)
        strandness = "+"
        outseqs = expandFusion(myTumorDNAsequence,  strandness ,  pepLen, leftLen)
        for aaseq in outseqs:
            if "*" not in aaseq:
                peptideIndex = "peptide_fusion_" + str(peptideID)
                print(peptideIndex, end="\t")
                tmpre = aaseq.split(":")   ## tmpre: peptide_155	LGLDWETQQ	aaCh	mutRegion	-	CTAGGGTTAGACTGGGAGACACAGCAG	17-43
                #print(tmpre[1] + "\t" + "aaCh" + "\t" + "mutRegion" + "\t" +  tmpre[2] + "\t" + tmpre[0] + "\t" + tmpre[3] + "\t" + tmpre[4] + "\t" + tmpre[5], end="\t")
                print("\t".join(tmpre), end="\t")
                print(">" + geneInfo)
                peptideID += 1
        #eventID += 1
        #subInfo = subSeq + ":" + subPep + ":" + rcStrandness + ":" + str(subSeqStart) + "-" + str(finalBase) + ":" + "reLoc=" + str(relativeLocRev) + ":" + eventLocation


