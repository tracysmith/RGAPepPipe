#!/usr/bin/env python
#This script will take a list of SRA filenames (w/o .sra extension) and output a VCF file for each.
#for each sample, a tab-delimited .txt file is required containing:
#SRAReadGroupID     SampleGroupID     LibraryID       SeqPlatform     paired/single
#ERR071082       ERS088906       ERX048850       illumina        paired    patient

import sys
import os
import argparse
from string import Template

    
def get_args():    
    parser = argparse.ArgumentParser(description='Pipeline for reference \
guided assembly')
    parser.add_argument("input", help="File describing read data information")
    parser.add_argument("reference", help="Fasta file of reference genome")
    parser.add_argument("dagtemplates", help="poolseq dag templates folder")

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
    return patDict

patDict = make_dict()

for pat in patDict:
    if len(patDict[pat]) > 1:
        print("There are {0} samples for patient {1}".format(len(patDict[pat]), pat))
        idList = patDict[pat]
        inputRealnList = []
        inputMpileupList = []
        realnTransferList = []
        mpileupTransferList = []
        indelsGtfList = []
        indelsGtfTransferList = []
        for i in idList:
            inputRealnList.append('-I {0}.realn.bam'.format(i))
            inputMpileupList.append('{0}.realn.patrealn.bam'.format(i))
            realnTransferList.append('{0}.realn.bam'.format(i))
            realnTransferList.append('{0}.realn.bai'.format(i))
            mpileupTransferList.append('{0}.realn.patrealn.bam'.format(i))
            mpileupTransferList.append('{0}.realn.patrealn.bai'.format(i))
            indelsGtfList.append('{0}_remIndels'.format(i))
            indelsGtfTransferList.append('{0}.indelreg.gtf'.format(i))
        variableMap = {}
        variableMap['ref'] = args.reference
        variableMap['pat'] = pat
        variableMap['inputs1'] = " ".join(inputRealnList)
        variableMap['trans1'] = ",".join(realnTransferList)
        variableMap['inputs2'] = " ".join(inputMpileupList)
        variableMap['trans2'] = ",".join(mpileupTransferList)
        variableMap['trans3'] = ",".join(indelsGtfTransferList)
        variableMap['indels'] = " ".join(indelsGtfList)
        print("Writing patient aware dag to {0}_pooled.dag".format(pat))
        with open("{0}_pooled.dag".format(variableMap['pat']), 'w') as dagfile:
            with open(args.dagtemplates + "patient_dag.template", 'r') as template_file:
                template = Template(template_file.read())
                out = template.substitute(variableMap)
                dagfile.write(out)
            for run in idList:
                dagfile.write('\n')
                variableMap['run'] = run
                with open(args.dagtemplates + "sample_pataware_dag.template", 'r') as template_file:
                    template = Template(template_file.read())
                    out = template.substitute(variableMap)
                    dagfile.write(out)
            dagfile.write("parent mpileup {0} child remMpileup".format(" ".join(indelsGtfList)))
    else:
        print("There is {0} samples for patient {1}".format(len(patDict[pat]), pat))
        run = patDict[pat][0]
        variableMap = {}
        variableMap['run'] = run
        variableMap['ref'] = args.reference
        print("Writing solo dag to {0}_pooled.dag".format(run))
        with open("{0}_pooled.dag".format(variableMap['run']), 'w') as dagfile:
            with open(args.dagtemplates + "sample_solo_dag.template", 'r') as template_file:
                template = Template(template_file.read())
                out = template.substitute(variableMap)
                dagfile.write(out)
