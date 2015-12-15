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
    parser.add_argument("-p", "--program", default="fastq_dump",
        help="program to start executing at: fastq_dump, fastqc, trim_galore, \
bwa (for longer reads), bwashort (for shorter reads), sort, dedup, readgroups, \
realign, vcf, cleanup [default: fastq_dump] (You must have run the pipeline up \
to the chosen starting program for this option to work)")
    parser.add_argument("-v", "--verbose", action="store_true", 
        help="Runs program in verbose mode which enables all debug output.")
    parser.add_argument("--ena", action="store_true", 
        help="Skips fastq_dump for every file.")
    parser.add_argument("-d","--fastqdir", default=".", 
        help="Specify a directory that contains .fastq files.")
    parser.add_argument("-t","--threads", type=int, default = 6,
        help="Specify the number of threads that should be used.")
    parser.add_argument("-s","--skip", action="store_true", 
        help="Continue even if errors occur.")

    return parser.parse_args()

args = get_args()

with open("rga_submit.template", "r") as template_file:
    template = Template(template_file.read())

with open(args.input, 'r') as infile:
    for line in infile:
        line = line.strip()
        inputList = line.split()
        variableMap = {}
        variableMap['ref'] = args.reference
        variableMap['run'] = inputList[0]
        variableMap['sample'] = inputList[1]
        variableMap['lib'] = inputList[2]
        variableMap['platform'] = inputList[3]
        pair = inputList[4]
        with open("{0}_rga.dag".format(variableMap['run']), "w") as dagfile:
            out = template.substitute(variableMap)
            dagfile.write(out)

    

