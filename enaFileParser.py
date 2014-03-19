#!/usr/bin/env python

import sys
import os

#This script takes a text file from the European Nucleotide Archive and produces
#two files. One is the correct input file for the SRAToVCF.py script. The other
#is a list of files to be downloaded from ENA.

# check for correct arguments
if len(sys.argv) != 2:
    print("Usage: enaFileParse.py <inputfile>")
    sys.exit(0)

inFile = open(sys.argv[1], 'r')
downloadFile = open(os.path.splitext(sys.argv[1])[0] + '_download.txt', 'w')
pipelineFile = open(os.path.splitext(sys.argv[1])[0] + '_pipelineIn.txt', 'w')

for i, line in enumerate(inFile):
    if i > 0:
        line = line.strip()
        entries = line.split('\t')
        readGroup = entries[5]
        sampleID = entries[3]
        experimentID = entries[4]
        platform = entries[7]
        if 'illumina' in platform.lower():
            platform = 'illumina'
        library = entries[8].lower()
        readURL = entries[9]
        pipelineFile.write('%s\t%s\t%s\t%s\t%s\n' % (readGroup, sampleID, 
                            experimentID, platform, library))
        for url in readURL.split(';'):
            downloadFile.write(url + '\n')
        
inFile.close()
downloadFile.close()
pipelineFile.close()
