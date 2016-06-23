#!/usr/bin/env python

import sys
import os
import argparse
from subprocess import call

# This script takes the realigned bam file from our pipeline and realigns it 
# with all other samples from the patient. It then produces an mpileup with
# options for analysis of pool-seq data.

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Realigns BAMs in a sample aware manner and gen pileups")
    input_file = parser.add_mutually_exclusive_group(required=True)
    input_file.add_argument('-i', '--inputFile',
        help ='pipeline infile with pat', 
        type=is_file)
    parser.add_argument('-p', '--pileup', action='store_true', 
        help='generate ind samp pileup')
    parser.add_argument('-m', '--mpileup', action='store_true',
        help='generate pat pileup')
    return parser.parse_args()

args = get_arguments()

def usage():
    print "funcEnr_makeFiles.py\n \
        -i <input file>\n \
        -p <generate ind pileup (default is None>\n \
        -m <generate mpileup (default is None)>"

def realign(inputList, pat):
    """Run realign for all samples of a patient"""
#    print("Processing sample " + pat)
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa -T RealignerTargetCreator -I ' + " -I ".join(inputList) + ' -o {pat}.intervals'.format(pat=pat), shell=True)
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa -T IndelRealigner -I ' + " -I ".join(inputList) + ' -targetIntervals {pat}.intervals -nWayOut .patrealn.bam'.format(pat=pat), shell=True)

def callable_loci(samp):
    """Run callable loci on samples"""
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I {samp}.realn.patrealn.bam -summary {samp}_defaults.summary -o {samp}_defaults.bed -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa'.format(samp=samp), shell=True)
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I {samp}.realn.patrealn.bam -summary {samp}_strict.summary -o {samp}_strict.bed -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa -frlmq 0.04 -mmq 20'.format(samp=samp), shell=True)

def ind_mpileup(samp):
    """Make mpileup for individual samples"""
    call('/opt/PepPrograms/RGAPipeline/samtools mpileup -B -q 20 -s -O -f /opt/data/mtuberculosis/MtbNCBIH37Rv.fa {samp}.realn.patrealn.bam > {samp}_Q20.pileup'.format(samp=samp), shell=True)
#    call('/opt/PepPrograms/RGAPipeline/samtools mpileup -B -u -v -S -q 20 -Q 20 -f /opt/data/mtuberculosis/MtbNCBIH37Rv.fa {samp}.realn.patrealn.bam > {samp}'.format(samp=samp), shell=True)

def rem_indels(samp):
    """Remove indels"""
    call('perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input {samp}_Q20.pileup --output {samp}_Q20.indelreg.gtf'.format(samp=samp), shell=True)
    call('perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {samp}_Q20.pileup --gtf {samp}_Q20.indelreg.gtf --output {samp}_Q20.noIndel.pileup'.format(samp=samp), shell=True)

def pat_mpileup(inputList, pat):
    """Make mpileup for patient"""
    print("Processing patient " + pat)
    print(" ".join(inputList))
#    call('/opt/PepPrograms/RGAPipeline/samtools mpileup -B -f /opt/data/mtuberculosis/MtbNCBIH37Rv.fa -q 20 -Q 20 ' + " ".join(inputList) + ' > {pat}.mpileup'.format(pat=pat), shell=True)
    call('/opt/PepPrograms/RGAPipeline/samtools mpileup -B -q 20 -f /opt/data/mtuberculosis/MtbNCBIH37Rv.fa ' + " ".join(inputList) + ' > {pat}.mpileup'.format(pat=pat), shell=True)

def rem_indels_mpileup(pat):
    """Remove indels"""
    call('perl /opt/PepPrograms/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input {pat}.mpileup --output {pat}.indelreg.gtf'.format(pat=pat), shell=True)
    call('perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {pat}.mpileup --gtf {pat}.indelreg.gtf --output {pat}.noIndel.mpileup'.format(pat=pat), shell=True)

def pat_sync(pat):
    """Make sync for patient"""
    print("Processing patient " + pat)
    call('java -ea -Xmx7g -jar /opt/PepPrograms/popoolation2_1201/mpileup2sync.jar --input {pat}.noIndel.mpileup --output {pat}_Q20.sync --fastq-type sanger --min-qual 20'.format(pat=pat), shell=True)
    #call('java -ea -Xmx7g -jar /home/peplab/src/popoolation2_1201/mpileup2sync.jar --input {pat}_filtered.mpileup --output {pat}_Q30.sync --fastq-type sanger --min-qual 30'.format(pat=pat), shell=True)

def get_RG(inFileName):
    patDict = {}
    with open(inFileName, 'r') as inFile:
        for line in inFile:
            RGID,RGSM,RGLB,RGPL,pair,patient=line.strip().split()
            if patient in patDict:
                patDict[patient].append(RGID)
            else:
                patDict[patient] = [RGID]
    print("There is/are " + str(len(patDict)) + " patient in this dataset.")
    for pat in patDict:
        print("There are " + str(len(patDict[pat])) + " samples for patient " + pat + ".") 
    return patDict


patDict = get_RG(args.inputFile)
for pat in patDict:
    idList = patDict[pat]
    inputRealnList = []
    inputMpileupList = []
    for i in idList:
        realn = i + '.realn.bam'
        patrealn = i + '.realn.patrealn.bam'
        inputRealnList.append(realn)
        inputMpileupList.append(patrealn)
    if args.mpileup is True:
        realign(inputRealnList, pat)
        pat_mpileup(inputMpileupList, pat)   
        rem_indels_mpileup(pat)
        pat_sync(pat)
    if args.pileup is True:
        map(ind_mpileup, idList)
        map(rem_indels, idList)
        map(callable_loci, idList)
