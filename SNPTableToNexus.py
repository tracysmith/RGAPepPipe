#!/usr/bin/env python

import sys
import re
import itertools

  
#SNP table from VCFstoSNPTable.py to nexus file 



#check for correct commandline arguments
if len(sys.argv) != 4 :
	print("Usage:  SNPTableToNexus.py  <input file>  <outputfile> \
               <name of reference sequence .fa file>")
	sys.exit(0)

InFileName = sys.argv[1]   #SNP table
OutFileName = sys.argv[2]  #Nexus format
Refname = sys.argv[3] #Reference sequence .fa
#NorALL = sys.argv[4] #Should all ambiguous bases be removed or just N

InFile = open(InFileName, 'r')	#access the file to be read from
OutFile = open(OutFileName, 'w') #overwrite what is present in the file
RefFile = open(Refname,'r')
PositionList = [] 
NameList = []
Refseq = ''
SeqList = []
linenumber = 0
SeqDict = {}

#generate reference sequence:
reflines = RefFile.readlines()
for line in reflines[1:len(reflines)]: #removes first line and combines lines into 1 string
    Refseq = Refseq+line.rstrip()

markerLine = '-'*len(Refseq)


#make a list of all the valid SNPs as lists, [position, SNP, sequence name], in the file
for Line in InFile:
    Line = Line.strip('\n')
    WordList = Line.split()
    SNP = WordList[3][3]
    position = int(WordList[2])-1
    name = WordList[0]

#fill dictionary
    if name not in SeqDict:
        SeqDict[name] = list(markerLine)
     
    SeqDict[name][position]=SNP #replace seq of "----" with "--S-" S->Any SNP



#write file header
OutFile.write("#NEXUS" + "\n" + "Begin DATA;" + "\n\t" + "Dimensions ntax=" +
	str(len(SeqDict)+1) + " nchar=" + str(len(Refseq)) + ";" + "\n\t" + "Format datatype=DNA gap=-;" 
	+ "\n\t" + "Matrix" + "\n")

#write ref sequence into file
OutFile.write(Refname + '\n')
OutFile.write("".join(Refseq))
OutFile.write('\n')

#write all of the SNPs into the new file
for key in SeqDict:
    newSeq = "".join(SeqDict[key])
    OutFile.write(key + "\n" + newSeq + "\n")
OutFile.write(";\nEND;\n")

InFile.close()
OutFile.close()
