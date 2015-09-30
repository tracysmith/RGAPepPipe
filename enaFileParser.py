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

rgIndex = 100
sampIndex = 100
exIndex = 100
platIndex = 100
libIndex = 100
urlIndex = 100

for i, line in enumerate(inFile):
    if i == 0:
        line = line.strip()
        entries = line.split('\t')
        for i,entry in enumerate(entries):
            if 'sample_accession' in entry:
                sampIndex = i
            elif 'experiment_accession' in entry:
                exIndex = i
            elif 'run_accession' in entry:
                rgIndex = i
            elif 'instrument' in entry:
                platIndex = i
            elif 'library_layout' in entry:
                libIndex = i
            elif 'fastq_ftp' in entry:
                urlIndex = i
            elif 'submitted_ftp' in entry:
                subIndex = i
        if 100 in [sampIndex, exIndex, rgIndex, platIndex, libIndex, urlIndex]:
            print "Cannot find all headers needed"
            sys.exit(2)
    else:
        line = line.strip()
        entries = line.split('\t')
        if len(entries) > urlIndex:
            readGroup = entries[rgIndex]
            sampleID = entries[sampIndex]
            experimentID = entries[exIndex]
            platform = entries[platIndex]
            if 'illumina' in platform.lower():
                platform = 'illumina'
            library = entries[libIndex].lower()
            readURL = entries[urlIndex]
            readSamp = entries[subIndex].split("/")[-1].split(".")[0]
            pipelineFile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (readGroup, 
            sampleID, experimentID, platform, library, readSamp))
            for url in readURL.split(';'):
                downloadFile.write(url + '\n')
        else:
            print('Problem with ' + entries[rgIndex])
inFile.close()
downloadFile.close()
pipelineFile.close()
