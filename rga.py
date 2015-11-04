#!/usr/bin/env python
#This script will take a list of SRA filenames (w/o .sra extension) and output a VCF file for each.
#for each sample, a tab-delimited .txt file is required containing:
#SRAReadGroupID     SampleGroupID     LibraryID       SeqPlatform     paired/single
#ERR071082       ERS088906       ERX048850       illumina        paired

import sys
import os
import subprocess
import shlex
import glob
import argparse

def call(cmd):
    #calls a system command with the subprocess module
    #redirects both stdout and stderr to a log file

    cmd = cmd.format(**(kvmap)) # replaces items in command with input values

    if args.verbose:
        print("Executing command: " + cmd + "\n")
    ret = subprocess.call(shlex.split(cmd))
    if(ret != 0):
        print("Pipeline did not complete successfully. \n Command : \n\n" + 
            cmd + "\n\n returned with non-zero code: " + str(ret))
        sys.exit(-1)
    
def call_bwa(cmd, filename):
    #bwa uses a redirect in its command so logging of stdout is unavailable
    cmd = cmd.format(**(kvmap))
    filename = filename.format(**(kvmap))
    
    bwa_out = open(filename, "w+")
    if args.verbose:
        print("Executing command: " + cmd + "\n")
    ret = subprocess.call(shlex.split(cmd), stdout=bwa_out)
    bwa_out.close()

def call_sort(cmd, cmd2):
    cmd = cmd.format(**(kvmap))
    cmd2 = cmd2.format(**(kvmap))
    
    if args.verbose:
        print("Executing command: " + cmd + "\n")
    samtools_view = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE) 
    subprocess.call(shlex.split(cmd2), stdin=samtools_view.stdout) 
    samtools_view.wait()

def fastq_dump():
    call("/opt/PepPrograms/sratoolkit.2.5.2-ubuntu64/bin/fastq-dump.2.5.2 \
--split-files {RGID}")
    
    if RGPL != "illumina":
        logfile = open(current_datetime+".log", "a+")
        logfile.write(line + '\n')
        logfile.write(RGID+ ' - not illumina' + '\n')
        logfile.close()
    if os.path.exists(fastqdir + '/' + RGID +'_2.fastq') != True and \
        pair == 'paired':
        call("rm {RGID}_1.fastq")
    else:
        fastqc()

def fastqc():
    filename = fastqdir + '/' + RGID + '_1.fastq'
    file = open(filename,'r')
    lines = file.readlines()
    fLine = lines[0]
    lposition = fLine.find('length')
    if lposition != -1:
        snum = fLine[lposition+7:len(fLine)-1]
    elif lposition == -1:
        snum = len(lines[1]) 
    numBP = int(snum)
    if pair == 'paired':
        call("/opt/PepPrograms/RGAPipeline/fastqc {fastqdir}/{RGID}_1.fastq \
{fastqdir}/{RGID}_2.fastq -t {threads} -o ./fastqc")
    elif pair == 'single':
        call("/opt/PepPrograms/RGAPipeline/fastqc {fastqdir}/{RGID}_1.fastq")
        
    trim_galore(numBP)

def trim_galore(numBP):
    if pair == 'paired':
        call("/opt/PepPrograms/RGAPipeline/trim_galore -q 15 --fastqc_args \
\"-t {threads} -o ./trimfastqc\"  -stringency 7 -o {pathToFastQ} --paired \
--retain_unpaired {fastqdir}/{RGID}_1.fastq {fastqdir}/{RGID}_2.fastq")
    elif pair == 'single':
        call("/opt/PepPrograms/RGAPipeline/trim_galore -q 15 --fastqc_args \
\"-t {threads} -o ./trimfastqc\"  -stringency 7 -o {pathToFastQ} {fastqdir}/\
{RGID}_1.fastq")
        
    if numBP < 70:
        bwa()
    else: #assumes numBP >= 70
        bwaMEM()

def bwaMEM():
    #Command to map with BWA (>70 bp)
    if pair =='paired':
        print("inside paired if")
        call_bwa("/opt/PepPrograms/RGAPipeline/bwa mem -M -t {threads} \
{reference} {pathToFastQ}/{RGID}_1_val_1.fq {pathToFastQ}/{RGID}_2_val_2.fq", 
            "{RGID}.sam")
    elif pair == 'single':
        call_bwa("/opt/PepPrograms/RGAPipeline/bwa mem -M -t {threads} \
{reference} {pathToFastQ}/{RGID}_1_trimmed.fq", "{RGID}.sam")
        
    sort()

def bwa():
    #Command to map with BWA (<70 bp)
    
    if pair == 'paired':
        call_bwa("/opt/PepPrograms/RGAPipeline/bwa aln -t {threads} \
{reference} {pathToFastQ}/{RGID}_1_val_1.fq", "{RGID}_1.sai")
        call_bwa("/opt/PepPrograms/RGAPipeline/bwa aln -t {threads} \
{reference} {pathToFastQ}/{RGID}_2_val_2.fq", "{RGID}_2.sai")
        call_bwa("/opt/PepPrograms/RGAPipeline/bwa sampe -P {reference} \
{RGID}_1.sai {RGID}_2.sai {pathToFastQ}/{RGID}_1_val_1.fq \
{pathToFastQ}/{RGID}_2_val_2.fq", "{RGID}.sam")
    elif pair == 'single':
        call_bwa("/opt/PepPrograms/RGAPipeline/bwa aln -t {threads} \
{reference} {pathToFastQ}/{RGID}_1_trimmed.fq", "{RGID}.sai")
        call_bwa("/opt/PepPrograms/RGAPipeline/bwa samse {reference} \
{RGID}.sai {pathToFastQ}/{RGID}_1_trimmed.fq", "{RGID}.sam")
        
    sort()

def sort():
    #Sort with samtools

    call_sort("/opt/PepPrograms/RGAPipeline/samtools view -bhSu \
{RGID}.sam", "/opt/PepPrograms/RGAPipeline/samtools sort - {RGID}.sort")
    
    dedup()
     
def dedup():
    #Command to mark duplicates in your bam file

    call("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/picard.jar \
MarkDuplicates I={RGID}.sort.bam O={RGID}.dedup.bam M={RGID}.metrics \
REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT")
    
    readgroups()
    
def readgroups():
    #Command to add/replace read groups

    call("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/picard.jar \
AddOrReplaceReadGroups I={RGID}.dedup.bam O={RGID}.ready.bam RGID={RGID} \
RGLB={RGLB} RGPL={RGPL} RGPU=dummy-barcode RGSM={RGSM} \
VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true")
    
    realign()

def realign():
    #Command to Realign with GATK
    
    call("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -I\
 {RGID}.ready.bam -R {reference} -T RealignerTargetCreator -o {RGID}.intervals") 
    call("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar \
-I {RGID}.ready.bam -R {reference} -T IndelRealigner -targetIntervals \
{RGID}.intervals -o {RGID}.realn.bam")
    
    vcf()
    
def vcf():
    #Command to Create VCF file

    call("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar \
-I {RGID}.realn.bam -R {reference} -T UnifiedGenotyper -o {RGID}.vcf \
-out_mode EMIT_ALL_CONFIDENT_SITES -stand_call_conf 20 -stand_emit_conf 20 \
--sample_ploidy 1 -nt {threads} -rf BadCigar")
    
    call("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysis\
TK.jar  -T VariantFiltration -R {reference} -o {RGID}_filter.vcf --variant \
{RGID}.vcf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum\
 < -12.5 || ReadPosRankSum < -8.0\" --filterName \"RGAPepPipeFilter\"")
    
    cleanup()

def cleanup():
    #Command to organize folder

    files = glob.glob(RGID + ".*");    
    files.extend(glob.glob(RGID + "_?.*"));
    
    call("mkdir -p {RGID}_output")

    for f in files:
        if f.endswith("_filter.vcf") or f.endswith(".realn.bam") or \
            f.endswith(".realn.bai"):
            call("mv {f} ./{RGID}_output".format(f = f))
        else:
            call("mv {f} ./intermediate_files".format(f = f))

    call("tar -czvf {RGID}.tar.gz {RGID}_output")
    
    
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

if args.ena:
    args.program = "fastqc"

program = args.program

inFileName = args.input
reference = args.reference
kvmap = {}

inFile = open(inFileName, 'r')

call("mkdir -p trim")
pathToFastQ = "trim"

fastqdir = os.path.abspath(args.fastqdir) 

call("mkdir -p trimfastqc")
call("mkdir -p fastqc")
call("mkdir -p intermediate_files")

#dictionary for starting the pipeline in the middle

functiondict =  {
        'fastq_dump':fastq_dump,
        'fastqc':fastqc,
        'trim_galore':trim_galore,
        'bwa':bwaMEM,
        'bwashort':bwa,
        'sort':sort,
        'dedup':dedup,
        'readgroups':readgroups,
        'realign':realign,
        'vcf':vcf,
        'cleanup':cleanup
    }

for line in inFile:
    numBP = 0
    line = line.strip()
    inputList = line.split()
    RGID = inputList[0]
    RGSM = inputList[1]
    RGLB = inputList[2]
    RGPL = inputList[3]
    pair = inputList[4]
    
    kvmap = {
        'RGID':RGID,
        'RGSM':RGSM,
        'RGLB':RGLB,
        'RGPL':RGPL,
        'reference':reference,
        'pathToFastQ':pathToFastQ,
        'fastqdir':fastqdir, 
        'threads':args.threads
    }
    
    #start the program at the specified program
    try:
        functiondict[program]()

    except IOError:
        print "\nFile {0} not found.\n".format(RGID)
        if(args.skip):
            continue
        else: 
            break
        
    except:
        if(args.skip):
            print "\nError occurred with file {0}, skipping.\n".format(RGID)
            continue
        else: 
            break


inFile.close()
