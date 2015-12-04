#!/usr/bin/env python

import sys
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
    """Create strings corresponding to contigs"""
    #contigs = []
    poscount = 0
    contig = ""
    secondLastLine = None
    lastLine = None
    with open(inFile, 'r') as vcf:
        for line in vcf:
            if ("#" in line):
                if ("CHROM" in line and "POS" in line and "ID" in line): 
                    line = line.strip()
                    tokens = line.split()
                    strain = tokens[9]
            else : #remaining lines
                line = line.strip()
                CHROM, POS, ID, REF, ALT, QUAL, FILTER = line.split()[0:7]
                if ALT == ".": #if not a SNP but confident site
                    ALLELE = REF
                elif ALT != "." and FILTER == "RGAPepPipeFilter":
                    ALLELE = "N"
                elif ALT != "." and FILTER == "PASS":
                    ALLELE = ALT
                if int(POS) == (poscount + 1):
                    contig = contig + ALLELE
                else:
                    contig = contig + "-"
                    #contigs.append(contig)
                    #contig = ALLELE
                poscount = int(POS) 
                secondLastLine = lastLine
                lastLine = line
    #contigs.append(contig)
    return(contig, strain) #contigs to contig

def write_fasta(contig, strain, RGID): #contigs to contig
    """Writes a dummy fasta alignment for each vcf"""
    outFile = strain + "_" + RGID + "_RGA.fasta"
    Sample = strain + "_" + RGID
    record = SeqRecord(Seq(contig, Gapped(IUPAC.ambiguous_dna, '-')), id=Sample)
    SeqIO.write(record, outFile, "fasta")
    #with open(outFile, 'w') as fasta:
    #    fasta.write(">" + Sample + '\n')
    #    fasta.write(contig + '\n')
        #for i, contig in enumerate(contigs):
        #    fasta.write(">contig" + str(i) + '\n')
        #    fasta.write(contig + '\n')


                    
for n in sys.argv[1:(len(sys.argv))] :
    RGID = n.split("_")[0]
    contig, strain = read_vcf(n) #contigs to contig
#    print("There are {0} contigs for {1}".format(len(contigs), RGID))
#    for i in contigs:
#        print len(i)
    write_fasta(contig, strain, RGID) #contigs to contig
