#!/usr/bin/env python
#This script will take a list of SRA filenames (w/o .sra extension) and output a VCF file for each.
#for each sample, a tab-delimited .txt file is required containing:
#SRAReadGroupID     SampleGroupID     LibraryID       SeqPlatform     paired/single
#ERR071082       ERS088906       ERX048850       illumina        paired

import sys
import os
import argparse
from string import Template

    
def get_args():    
    parser = argparse.ArgumentParser(description='Pipeline for reference \
guided assembly')
    parser.add_argument("input", help="File describing read data information")
    parser.add_argument("reference", help="Fasta file of reference genome")
    parser.add_argument("dagtemplate", help="dag template")

    return parser.parse_args()

args = get_args()

def make_dict():
    patDict = {}
    with open(args.input, 'r') as infile:
        for line in infile:
            line = line.strip().split()
            RGID = line[0]
            RGSM = line[1]
            strat = line[4]
            pat = line[5]
            if pat in patDict:
                patDict[pat].append(RGID)
            else:
                patDict[pat] = [RGID]
        print("There is/are {0} patient(s) in this dataset.".format(len(patDict)))
        for pat in patDict:
            print("There is/are {0} samples for patient {1}.".format(len(patDict[pat]), pat))
    return patDict

patDict = make_dict()

with open(args.dagtemplate, "r") as template_file:
    template = Template(template_file.read())

for pat in patDict:
    idList = patDict[pat]
    inputRealnList = []
    inputMpileupList = []
    realnTransferList = []
    mpileupTransferList = []
    for i in idList:
        inputRealnList.append('-I {0}.realn.bam'.format(i))
        inputMpileupList.append('{0}.realn.patrealn.bam'.format(i))
        realnTransferList.append('{0}.realn.bam'.format(i))
        realnTransferList.append('{0}.realn.bai'.format(i))
        mpileupTransferList.append('{0}.realn.patrealn.bam'.format(i))
        mpileupTransferList.append('{0}.realn.patrealn.bai'.format(i))
    variableMap = {}
    variableMap['ref'] = args.reference
    variableMap['pat'] = pat
    variableMap['inputs1'] = " ".join(inputRealnList)
    variableMap['trans1'] = ",".join(realnTransferList)
    variableMap['inputs2'] = " ".join(inputMpileupList)
    variableMap['trans2'] = ",".join(mpileupTransferList)
    with open("{0}_pooled.dag".format(variableMap['pat']), 'w') as dagfile:
        out = template.substitute(variableMap)
        dagfile.write(out)
        for run in idList:
            dagfile.write("\nJOB {0}_pileup pileup.submit\n".format(run))
            dagfile.write('VARS {0}_pileup REF="{1}"\n'.format(run, args.reference))
            dagfile.write('VARS {0}_pileup RUN="{0}"\n\n'.format(run))
            dagfile.write('JOB {0}_remIndels remIndels.submit\n'.format(run))
            dagfile.write('VARS {0}_remIndels RUN="{0}"\n\n'.format(run))
            dagfile.write('parent realign child {0}_pileup\n'.format(run))
            dagfile.write('parent {0}_pileup child {0}_remIndels\n'.format(run))
