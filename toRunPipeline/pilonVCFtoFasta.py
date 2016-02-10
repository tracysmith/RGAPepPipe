#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped


#####################
# This script takes a variable number of input vcf files.  Each vcf file must contain
# snp information about only one strain.  At least one input file must be supplied.  The
# name of an output file must also be supplied, and output is in snpt able format.  The 
# output contains the combined information about all strains.
# NOTE : position for indels is the locus before the insertion or deletion
#############


# check for correct number of inputs
if len(sys.argv) < 2 :
     print("Usage: VCFToFasta.py <input vcf 1> ... <input vcf n>")
     sys.exit(0)
    
def read_vcf(inFile):
    """Create strings corresponding to chromosome"""
    with open(inFile, 'r') as vcf:
        for line in vcf:
            if ("contig" in line and "length" in line):
                line = line.strip()
                length = line.split("=")[-1].strip(">")
                print length
                refID = line.split(",")[0].split("=")[-1]
                chromosome = ["?"] * int(length)
            #elif ("CHROM" in line and "POS" in line and "ID" in line): 
            #    line = line.strip()
            #    tokens = line.split()
            #    strain = tokens[9]
            elif line[0] != "#":
                line = line.strip()
                CHROM, POS, ID, REF, ALT, QUAL, FILTER = line.split()[0:7]
                FILTERS = FILTER.split(";")
                if len(REF) == 1 and len(ALT) == 1:
                    if "Amb" in FILTERS:
                        ALLELE = "?"
                    elif "PASS" in FILTERS:
                        if ALT == ".":
                            ALLELE = REF
                        else:
                            ALLELE = ALT
                    elif "Del" in FILTERS:
                            ALLELE = "-"
                    else:
                        if FILTER == "LowCov":
                            continue
                        else:
                            print(FILTER)
                index = int(POS) - 1
                chromosome[index] = ALLELE
    return(chromosome, refID)

def write_fasta(chromosome, RGID, refID):
    """Writes a RGA fasta alignment for each vcf"""
    outFile = RGID + "_RGA_pilon.fasta"
    Sample = "pilon_" + RGID
    record = SeqRecord(Seq("".join(chromosome), Gapped(IUPAC.ambiguous_dna, '-')), id=Sample, description = "RGA_to_" + refID)
    SeqIO.write(record, outFile, "fasta")

for n in sys.argv[1:(len(sys.argv))] :
    RGID = os.path.basename(n).split("_")[0]
    chromosome, refID = read_vcf(n)
    write_fasta(chromosome, RGID, refID)
