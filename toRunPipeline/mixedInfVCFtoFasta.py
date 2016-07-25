#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio import AlignIO
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

def find_regions():
    regions = []
    with open("/opt/data/mtuberculosis/160129_Mtb_removeRegions_mergedCloserThan1000bp.bed", 'r') as remReg:
        for line in remReg:
            line = line.strip().split('\t')
            start = int(line[1])
            stop = int(line[2])
            regions.append((start,stop))
    return regions

def read_vcf(inFile):
    """Create strings corresponding to chromosome"""
    d = {}
    #outFile = RGID + "_mixed.txt"
    with open(inFile, 'r') as vcf: #, open(outFile, 'w') as out:
        #out.write("pos\tA\tC\tG\tT\n")
        for line in vcf:
            if ("contig" in line and "length" in line):
                line = line.strip()
                length = line.split("=")[-1].strip(">")
                print length
                refID = line.split(",")[0].split("=")[-1]
                #geno1 = ["?"] * int(length)
                #geno2 = ["?"] * int(length)
            #elif ("CHROM" in line and "POS" in line and "ID" in line): 
            #    line = line.strip()
            #    tokens = line.split()
            #    strain = tokens[9]
            elif line[0] != "#":
                line = line.strip()
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = line.split()[0:8]
                FILTERS = FILTER.split(";")
                if len(REF) == 1 and len(ALT) == 1:
                    if "Del" in FILTERS:
                        ALLELE = "NA,NA,NA,NA"
                    elif "LowCov" in FILTERS:
                        ALLELE = "NA,NA,NA,NA"
                    else:
                        fields = INFO.split(";")
                        for i in fields:
                            if "QP" in i:
                                ALLELE = i.split("=")[1]
                    alleles = ALLELE.split(",")
                d[int(POS)] = alleles
                #out.write(POS + '\t' + "\t".join(alleles) + '\n')
                #index = int(POS) - 1
                #chromosome[index] = ALLELE
    #return(chromosome, refID)
    return d

def remove_regions(d, regions):
    """Convert PPE/rep regions to NAs"""
    for pos in d:
        for reg in regions:
            if int(pos) >= reg[0] and int(pos) <= reg[1]:
                d[pos] = ["NA", "NA", "NA", "NA"]
                #print("PPE REGION")
            else:
                continue
    return d

def write_txt(d):
    """Write text file with counts for each position"""
    outFile = RGID + "_mixed.txt"
    with open(outFile, 'w') as out:
        out.write("pos\tA\tC\tG\tT\n")
        for i in range(1,4411533):
            if i in d:
                out.write(str(i) + '\t' + "\t".join(d[i]) + '\n')

def write_fasta():
    """Writes a RGA fasta alignment for each vcf"""
    outFile = RGID + "_mixed.fasta"
    geno1 = ["-"] * 4411532 
    geno2 = ["-"] * 4411532
    gt = {}
    gt[0] = 'A'
    gt[1] = 'C'
    gt[2] = 'G'
    gt[3] = 'T'
    with open(RGID + '_mixed.txt', 'r') as inFile:
        next(inFile)
        for line in inFile:
            line = line.strip().split('\t')
            pos = int(line[0])
            print pos
            try:
                counts = [int(i) for i in line[1:4]]
                majAF = float(max(counts))/sum(counts)
                print majAF
                majAL = counts.index(max(counts))
                print gt[majAL]
                minAL = counts.index(sorted(counts, reverse=True)[1])
                print gt[minAL]
                if majAF >= 0.55 and majAF <= 0.65:
                    print("mixed genotype")
                    idx = pos - 1
                    geno1[idx] = gt[majAL]
                    geno2[idx] = gt[minAL]
                elif majAF >= 0.95:
                    idx = pos - 1
                    geno1[idx] = gt[majAL]
                    geno2[idx] = gt[majAL]
                else:
                    continue
            except ZeroDivisionError:
                continue
            except ValueError:
                continue
    record1 = SeqRecord(Seq("".join(geno1), Gapped(IUPAC.ambiguous_dna, '-')), id=RGID + "_Geno1", description = "mixed infection")
    record2 = SeqRecord(Seq("".join(geno1), Gapped(IUPAC.ambiguous_dna, '-')), id=RGID + "_Geno2", description = "mixed infection")
    genotypes = []
    genotypes.append(record1)
    genotypes.append(record2)
    SeqIO.write(genotypes, outFile, "fasta")
    return geno1, geno2

for n in sys.argv[1:(len(sys.argv))] :
    RGID = os.path.basename(n)[0:-4]
    #d = read_vcf(n)
    #regions = find_regions()
    #d = remove_regions(d, regions)
    #write_txt(d)
    geno1, geno2 = write_fasta()
    #chromosome, refID = read_vcf(n)
    #write_fasta(chromosome, RGID, refID)
