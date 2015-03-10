#!/usr/bin/python
#; -*- mode: python -*- #uncomment this when ready to run on server
#This script will take a list of SRA filenames (w/o .sra extension) and output a VCF file for each.
#for each sample, a tab-delimited .txt file is required containing:
#SRAReadGroupID     SampleGroupID     LibraryID       SeqPlatform     paired/single
#ERR071082       ERS088906       ERX048850       illumina        paired
#Must be run from folder with indexed reference genome files in it as well as downloaded .sra files
#When running the script, type:
#python [path to script file] [path to .txt file] [.fasta reference file] [options]


import sys,os,subprocess, shlex, glob
from datetime import datetime
from optparse import OptionParser

current_datetime = datetime.today().strftime("%d-%m-%Y-%H%M")

#calls a system command with the subprocess module
#redirects both stdout and stderr to a log file
def call_with_log(cmd):
    cmd = cmd.format(**(kvmap))

    logfile = open(current_datetime+".log", "a+")
    if options.verbose:
        print("Executing command: " + cmd + "\n")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=logfile, stderr=logfile)
    if(ret != 0):
        print("Pipeline did not complete successfully. \n Command : \n\n" + cmd + "\n\n returned with non-zero code: " + str(ret))
        logfile.write("Pipeline did not complete successfully. \n Command : \n\n" + cmd + "\n\n returned with non-zero code: " + str(ret))
        logfile.close()
        sys.exit(-1)
    logfile.close()
    
#bwa uses a redirect in its command so logging of stdout is unavailable 
def call_bwa_with_log(cmd, filename):
    cmd = cmd.format(**(kvmap))
    filename = filename.format(**(kvmap))
    
    bwa_out = open(filename, "w+")
    logfile = open(current_datetime+".log", "a+")
    if options.verbose:
        print("Executing command: " + cmd + "\n")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=bwa_out, stderr=logfile)
    logfile.close()
    bwa_out.close()

#sort uses a pipe so it must also be modified
def call_sort_with_log(cmd, cmd2):
    cmd = cmd.format(**(kvmap))
    cmd2 = cmd2.format(**(kvmap))
    
    logfile = open(current_datetime+".log", "a+")
    if options.verbose:
        print("Executing command: " + cmd + "\n")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    samtools_view = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=logfile)
    subprocess.call(shlex.split(cmd2), stdin=samtools_view.stdout, stderr=logfile)
    samtools_view.wait()
    logfile.close()

def fastq_dump():
    print("fastq_dump started")
    
    call_with_log("/opt/PepPrograms/RGAPipeline/fastq-dump --split-files {RGID}.sra")
    
    print("fastq_dump completed")
    
    if RGPL != "illumina":
        logfile = open(current_datetime+".log", "a+")
        logfile.write(line + '\n')
        logfile.write(RGID+ ' - not illumina' + '\n')
        logfile.close()
    if os.path.exists(fastqdir + '/' + RGID +'_2.fastq') != True and pair == 'paired':
        call_with_log("rm {RGID}_1.fastq")
    else:
        fastqc()

def fastqc():
    print("fastqc started")
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
        call_with_log("/opt/PepPrograms/RGAPipeline/fastqc {fastqdir}/{RGID}_1.fastq {fastqdir}/{RGID}_2.fastq -t {threads} -o ./fastqc")
    elif pair == 'single':
        call_with_log("/opt/PepPrograms/RGAPipeline/fastqc {fastqdir}/{RGID}_1.fastq")
        
    print("fastqc completed")
    
    trim_galore(numBP)

def trim_galore(numBP):
    print("trim_galore started")
    
    if pair == 'paired':
        call_with_log("/opt/PepPrograms/RGAPipeline/trim_galore -q 15 --fastqc_args \"-t {threads} -o ./trimfastqc\"  -stringency 7 -o {pathToFastQ} --paired --retain_unpaired {fastqdir}/{RGID}_1.fastq {fastqdir}/{RGID}_2.fastq")
    elif pair == 'single':
        call_with_log("/opt/PepPrograms/RGAPipeline/trim_galore -q 15 --fastqc_args \"-t {threads} -o ./trimfastqc\"  -stringency 7 -o {pathToFastQ}  {fastqdir}/{RGID}_1.fastq")
        
    print('trim_galore completed')
    
    if numBP >= 70:
        bwaMEM()
    else:
        bwa()

#Command to map with BWA (>70 bp)
def bwaMEM():
    print("bwa started")
    
    if pair =='paired':
        print("inside paired if")
        call_bwa_with_log("/opt/PepPrograms/RGAPipeline/bwa mem -M -t {threads} {reference} {pathToFastQ}/{RGID}_1_val_1.fq {pathToFastQ}/{RGID}_2_val_2.fq", "{RGID}.sam")
    elif pair == 'single':
        call_bwa_with_log("/opt/PepPrograms/RGAPipeline/bwa mem -M -t {threads} {reference} {pathToFastQ}/{RGID}_1_trimmed.fq", "{RGID}.sam")
        
    print('bwa completed')
    
    if os.path.exists(RGID+'.sam') == True:
        sort()

#Command to map with BWA (<70 bp)
def bwa():
    print("bwa started")
    
    if pair == 'paired':
        call_bwa_with_log("/opt/PepPrograms/RGAPipeline/bwa aln -t {threads} {reference} {pathToFastQ}/{RGID}_1_val_1.fq", "{RGID}_1.sai")
        call_bwa_with_log("/opt/PepPrograms/RGAPipeline/bwa aln -t {threads} {reference} {pathToFastQ}/{RGID}_2_val_2.fq", "{RGID}_2.sai")
        call_bwa_with_log("/opt/PepPrograms/RGAPipeline/bwa sampe -P {reference} {RGID}_1.sai {RGID}_2.sai {pathToFastQ}/{RGID}_1_val_1.fq {pathToFastQ}/{RGID}_2_val_2.fq", "{RGID}.sam")
    elif pair == 'single':
        call_bwa_with_log("/opt/PepPrograms/RGAPipeline/bwa aln -t {threads} {reference} {pathToFastQ}/{RGID}_1_trimmed.fq", "{RGID}.sai")
        call_bwa_with_log("/opt/PepPrograms/RGAPipeline/bwa samse {reference} {RGID}.sai {pathToFastQ}/{RGID}_1_trimmed.fq", "{RGID}.sam")
        
    print('bwa completed')
    
    if os.path.exists(RGID+'.sam') == True:
        sort()

#Sort with samtools
def sort():
    print("sort started")
    
    call_sort_with_log("/opt/PepPrograms/RGAPipeline/samtools view -bhSu {RGID}.sam", "samtools sort - {RGID}.sort")
    
    print('sort completed')
    
    dedup()
     
#Command to mark duplicates in your bam file
def dedup():
    print("dedup started")
    
    call_with_log("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/MarkDuplicates.jar I={RGID}.sort.bam O={RGID}.dedup.bam M={RGID}.metrics REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT")
    
    print('dedup completed')
    
    readgroups()
    
#Command to add/replace read groups
def readgroups():
    print("readgroups started")
    
    call_with_log("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/AddOrReplaceReadGroups.jar I={RGID}.dedup.bam O={RGID}.ready.bam RGID={RGID} RGLB={RGLB} RGPL={RGPL} RGPU=dummy-barcode RGSM={RGSM} VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true")
    
    print('readgroups completed')
    
    realign()

#Command to Realign with GATK
def realign():
    print("realign started")
    
    call_with_log("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -I {RGID}.ready.bam -R {reference} -T RealignerTargetCreator -o {RGID}.intervals") 

    call_with_log("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -I {RGID}.ready.bam -R {reference} -T IndelRealigner -targetIntervals {RGID}.intervals -o {RGID}.realn.bam")
    
    print('realign completed')
    
    vcf()
    
#Command to Create VCF file
def vcf():
    print("vcf started")
    
    call_with_log("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -I {RGID}.realn.bam -R {reference} -T UnifiedGenotyper -o {RGID}.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -stand_call_conf 20 -stand_emit_conf 20 --sample_ploidy 1 -nt {threads} -rf BadCigar")
    call_with_log("java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -I {RGID}.realn.bam -R {reference} -T DepthOfCoverage -L {RGID}.intervals -U -S SILENT -rf BadCigar")
    
    print('vcf completed')
    
    if os.path.exists (RGID+'.vcf') == True:
        cleanup()

#Command to organize folder
def cleanup():
    print("cleanup started")
    
    call_with_log("mv {RGID}.vcf {projectname}_vcf")
    
    files = glob.glob(RGID + ".*");    
    files.extend(glob.glob(RGID + "_?.*"));
    
    for file in files:
        if file.endswith(".fastq") or file.endswith(".realn.bam") or file.endswith(".realn.bai"):
            call_with_log("mv {file} ./saved_files".format(file = file) )

        else:
            call_with_log("mv {file} ./intermediate_files".format(file = file))
   
    print ('cleanup completed')

    
##begin main execution


##for command line argument processing
    
parser = OptionParser(usage="usage: %prog [options] <path to .txt file> <.fasta reference file>")

parser.add_option("-p", "--program", dest="start_prog", help="program to start executing at: fastq_dump, fastqc, trim_galore, bwa, sort, dedup, readgroups, realign, vcf, cleanup [default: fastq_dump] (You must have run the pipeline up to the chosen starting program for this option to work)")

parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Runs program in verbose mode which enables all debug output.")

parser.add_option("--ena", action="store_true", dest="ena", help="Skips fastq_dump for every file. For read data already in fastq format.")

parser.add_option("-d","--fastqdir", dest="fastqdir",default=".", help="Specify a directory that contains .fastq files.")

parser.add_option("-t","--threads", type="int", dest="numthreads", help="Specify the number of threads that should be used.")

parser.set_defaults(start_prog="fastq_dump", verbose=False, ena=False, numthreads=6)

(options, args) = parser.parse_args()

if(len(args) < 2):
    print("usage: <path to .txt file> <.fasta reference file> [options]")
    sys.exit(-1)

if options.ena:
    options.start_prog="fastqc"

program = options.start_prog

inFileName = args[0]
reference = args[1]


inFile = open(inFileName, 'r')
projectname = os.path.basename( inFileName.strip(".txt") )
kvmap= {'projectname':projectname}

if(len(args) >= 3):
    pathToFastQ = os.path.abspath(args[-1])
    if not os.path.isdir(pathToFastQ):
        call_with_log("mkdir -p " + pathToFastQ)
else:
    call_with_log("mkdir -p trim")
    pathToFastQ = "trim"

fastqdir = os.path.abspath(options.fastqdir) #Set fastqdir to what was specified, or "" if not

call_with_log("mkdir -p trimfastqc")
call_with_log("mkdir -p fastqc")
call_with_log("mkdir -p {projectname}_vcf")
call_with_log("mkdir -p intermediate_files")
call_with_log("mkdir -p saved_files")

#dictionary for starting the pipeline in the middle

functiondict =  {
        'fastq_dump':fastq_dump,
        'fastqc':fastqc,
        'trim_galore':trim_galore,
        'bwa':bwa,
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
        'projectname':projectname,
        'fastqdir':fastqdir, 
        'threads':options.numthreads
    }
    
    #start the program at the specified program
    functiondict[options.start_prog]()
    if options.ena:
        options.start_prog="fastqc"
    else:
        #options.start_prog = "fastq_dump"
        options.start_prog = program
            

inFile.close()

#following code is for debugging purposes only:        
"""for line in InFile:
    line = line.strip()
    inputList = line.split()
    RGID = inputList[0]
    RGSM = inputList[1]
    RGLB = inputList[2]
    RGPL = inputList[3]
    pair = inputList[4]

    bwa()
    if os.path.exists(RGID+'.sam') == True:
        print('bwa completed')
        sort()
        print('sort completed')
        dedup()
        print('dedup completed')
        readgroups()
        print('readgroups completed')
        realign()
        print('realign completed')
        vcf()
        print('vcf completed')
        if os.path.exists (RGID+'.vcf') == True:
            print('vcf completed')
            cleanup()
            print ('cleanup completed')
        else:
            failedF.write(RGID + ' - vcf' + '\n')
    else:
        failedF.write(RGID + ' - bwa' + '\n')"""

