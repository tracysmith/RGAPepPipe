#!/bin/usr/env python

import sys

#####################
# This script takes a variable number of input vcf files.  Each vcf file must contain
# snp information about only one strain.  At least one input file must be supplied.  The
# name of an output file must also be supplied, and output is in snpt able format.  The 
# output contains the combined information about all strains.
# NOTE : position for indels is the locus before the insertion or deletion
#############

snpDict = {}
locs = []

# check for correct number of inputs
if len(sys.argv) < 3 :
     print("Usage: VCFsToSnpTable.py <input vcf 1> ... <input vcf n> <output file>")
     sys.exit(0)

outFile = open(sys.argv[(len(sys.argv) - 1)], 'w')
for n in sys.argv[1:(len(sys.argv) - 1)] :
    
    inFile = open(n, "r")

    # parse data into a dictionary indexed by locus of snp

    lines = inFile.readlines()
    print(len(lines))

    for line in lines :
        if ("#" in line) :
            if ("CHROM" in line and "POS" in line and "ID" in line) : #reads first line of values to get strain
                line = line.strip()
                tokens = line.split()
                strain = tokens[9]
        else : #remaining lines
            line = line.strip()
            tokens = line.split()
            # each tuple contains strain name, reference base and snp
            if tokens[4] == ".": #if not a SNP but confident site
                tokens[4] = tokens[3]
            tup = strain, tokens[3], tokens[4]
            outFile.write(tup[0] + "\tlocus " + tokens[1] + "\t" + tup[1] + "->" + tup[2] + "\n")

    inFile.close()

outFile.close()

