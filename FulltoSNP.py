#!/usr/bin/env python

import sys
import re
import itertools
import math
from Bio import SeqIO

#SNP alignment from full Alignment nexus file

#Check for correct commandline arguments
if len(sys.argv) != 4:
    print("Usage: FulltoSNP.py <nexus file> <output file> <threshold>")
    sys.exit(0)

#Get filenames
InFileName = sys.argv[1]
OutFileName = sys.argv[2]
threshold = sys.argv[3]
PosOutFileName = sys.argv[2]+'positions'

record_dict = SeqIO.to_dict(SeqIO.parse(InFileName,"fasta"))

#seperate speciesnames from sequences
seqs = []
titles = []

for key in record_dict:
    titles.append(key)
    x = record_dict[key]
    seqs.append(x.seq)



#transpose string lists
thresh = math.ceil(float(threshold) * len(seqs))
print(thresh)
seqsTran = zip(*seqs)
snps = []
#for every tuple check if value is the same, if so remove tuple
pos = 1
positions=[]
for s in seqsTran[:]:
    if len(set(s))!=1 and s.count('-')<= thresh:
	snps.append(s)
	positions.append(pos)
    pos=pos+1
print(len(positions))
seqsTran = []
results = zip(*snps)
for i in range(len(results)):
	results[i] = ''.join(results[i])

SeqDict={}
print(len(results[0]))
for i in range(len(results)):
    SeqDict[titles[i]]=results[i]

OutFile = open(OutFileName,'w')

#write file header
OutFile.write("#NEXUS" + "\n" + "Begin DATA;" + "\n\t" + "Dimensions ntax=" + str(len(SeqDict)) + " nchar=" + str(len(results[0])) + ";" + "\n\t" + "Format datatype=DNA gap=-;" + "\n\t" + "Matrix" + "\n")



#write all of the SNPs into the new file
for key in SeqDict:
    newSeq = "".join(SeqDict[key])
    OutFile.write("'" + key + "'" + "\t" + newSeq + "\n")
OutFile.write(";" + "\n" + "END;")
OutFile.close()
OutFile2 = open(PosOutFileName,'w')
for i in positions:
	OutFile2.write(str(i)+'\n')
OutFile2.close()
