#!/usr/bin/python
#; -*- mode: python -*- #uncomment this when ready to run on server
#This script will take a list of SRA filenames (w/o .sra extension) and output a VCF file for each.
#for each sample, a tab-delimited .txt file is required containing:
#SRAReadGroupID	 SampleGroupID	 LibraryID       SeqPlatform	 paired/single
#ERR071082       ERS088906       ERX048850       illumina        paired
#Must be run from folder with indexed reference genome files in it as well as downloaded .sra files
#Must also have a folder within this folder called "trimfastqc"
#When running the script, type:
#ipython [path to script file] [path to .txt file] [.fasta reference file] [path to directory named "trim"]


import sys,os,subprocess
from datetime import datetime

current_datetime = datetime.today().strftime("%d-%m-%Y-%H%M")

#calls a system command with the subprocess module
#redirects both stdout and stderr to a log file
def call_with_log(cmd):
    logfile = open(current_datetime+".log", "a+")
    subprocess.call(cmd, stdout=logfile, stderr=logfile)
    logfile.close()

def fastq_dump():
    call_with_log(["/opt/PepPrograms/RGAPipeline/fastq-dump", "--split-files" + RGID + ".sra"])

def fastqc():
    filename = RGID+ '_1.fastq'
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
        call_with_log(["/opt/PepPrograms/RGAPipeline/fastqc", RGID+"_1.fastq", RGID+"_2.fastq", "-t 6", "-o ./fastqc"])
    elif pair == 'single':
        call_with_log(["/opt/PepPrograms/RGAPipeline/fastqc", RGID+"_1.fastq"])
    return numBP

def trim_galore():
    if pair == 'paired':
        call_with_log(["/opt/PepPrograms/RGAPipeline/trim_galore", "-q 15", "--fastqc_args \"-t 6 -o ./trimfastqc\"", "-stringency 7", "-o ./trim", "--paired", "--retain_unpaired", RGID+"_1.fastq", RGID+"_2.fastq"])
    elif pair == 'single':
        call_with_log(["/opt/PepPrograms/RGAPipeline/trim_galore", "-q 15", "--fastqc_args \"-t 6 -o ./trimfastqc\"",  "-stringency 7", "-o ./trim",  RGID+"_1.fastq"])
      

#Command to map with BWA (>70 bp)
def bwaMEM():
    print("inside mem")
    if pair =='paired':
	print("inside paired if")
        call_with_log(["/opt/PepPrograms/RGAPipeline/bwa", "mem", "-M", "-t 6", reference, pathToFastQ + "/" + RGID + "_1_val_1.fq", "pathToFastQ/"+ RGID + "_2_val_2.fq", "> " + RGID + ".sam"])
    elif pair == 'single':
        call_with_log(["/opt/PepPrograms/RGAPipeline/bwa", "mem", "-M", "-t 6", reference, pathToFastQ + "/" + RGID + "_1_trimmed.fq",  "> " + RGID + ".sam"])

#Command to map with BWA (<70 bp)
def bwa():
    if pair == 'paired':
        call_with_log(["/opt/PepPrograms/RGAPipeline/bwa", "aln", "-t 6", reference,  pathToFastQ + "/" + RGID + "_1_val_1.fq", "> " + RGID + "_1.sai"])
        call_with_log(["/opt/PepPrograms/RGAPipeline/bwa", "aln", "-t 6", reference,  pathToFastQ + "/" + RGID + "_2_val_2.fq", "> " + RGID + "_2.sai"])
        call_with_log(["/opt/PepPrograms/RGAPipeline/bwa", "sampe", "-P", reference,  RGID +  "_1.sai", RGID + "_2.sai", pathToFastQ + "/" + RGID + "_1_val_1.fq", pathToFastQ + "/" + RGID + "_2_val_2.fq", "> " + RGID + ".sam"])
    elif pair == 'single':
        call_with_log(["/opt/PepPrograms/RGAPipeline/bwa", "aln", "-t 6", reference, pathToFastQ + "/" + RGID + "_1_trimmed.fq", "> " + RGID + ".sai"])
        call_with_log(["/opt/PepPrograms/RGAPipeline/bwa", "samse", reference, RGID + ".sai", pathToFastQ + "/" + RGID + "_1_trimmed.fq", "> " + RGID + ".sam"])

#Sort with samtools
def sort():
    call_with_log(["/opt/PepPrograms/RGAPipeline/samtools", "view", "-bhSu", RGID + ".sam", "|", "samtools", "sort", "- " + RGID + ".sort"])

#Command to mark duplicates in your bam file
def dedup():
    call_with_log(["java", "-Xmx2g", "-jar /opt/PepPrograms/RGAPipeline/MarkDuplicates.jar", "I=" + RGID + ".sort.bam", "O=" + RGID + ".dedup.bam", "M=" + RGID + ".metrics", "REMOVE_DUPLICATES=true", "AS=true", "VALIDATION_STRINGENCY=SILENT"])

#Command to add/replace read groups
def readgroups():
    call_with_log(["java", "-Xmx2g", "-jar /opt/PepPrograms/RGAPipeline/AddOrReplaceReadGroups.jar", "I=" + RGID + ".dedup.bam", "O=" + RGID + ".ready.bam" "RGID=" + RGID, "RGLB=" + RGLB, "RGPL=" + RGPL, "RGPU=dummy-barcode", "RGSM=" + RGSM, "VALIDATION_STRINGENCY=SILENT", "SORT_ORDER=coordinate", "CREATE_INDEX=true"])

#Command to Realign with GATK
def realign():
    call_with_log(["java", "-jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar", "-I " +  RGID + ".ready.bam", "-R " + reference, "-T RealignerTargetCreator", "-o "  + RGID + ".intervals"]) 

    call_with_log(["java", "-jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar", "-I " +  RGID + ".ready.bam", "-R " + reference, "-T IndelRealigner", "-targetIntervals " + RGID + ".intervals", "-o " + RGID + ".realn.bam"])
	
#Command to Create VCF file
def vcf():
    call_with_log(["java", "-jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar", "-I " + RGID + ".realn.bam", "-R " + reference, "-T UnifiedGenotyper", "-o " + RGID + ".vcf", "-out_mode EMIT_ALL_CONFIDENT_SITES", "-stand_call_conf 20", "-stand_emit_conf 20", "--sample_ploidy 1", "-nt 6", "-rf BadCigar"])
    call_with_log(["java", "-jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar", "-I " + RGID + ".realn.bam", "-R " + reference, "-T DepthOfCoverage", "-L " + RGID + ".intervals", "-U", "-S SILENT", "-rf BadCigar"])

#Command to organize folder
def cleanup():
    call_with_log(["mv", RGID + ".vcf", "./" + projectname + "_vcf"])
    call_with_log(["mv", RGID + ".*", "./intermediate_files"])
    call_with_log(["mv", RGID + "_?.*",  "./intermediate_files"])

if(len(sys.argv) != 4):
    print("usage: ipython <path to script file> <path to .txt file> <.fasta reference file> <path to directory named 'trim'>")
    sys.exit(-1)

inFileName = sys.argv[1]
reference = sys.argv[2]
pathToFastQ = sys.argv[3]

inFile = open(InFileName, 'r')

projectname = InFileName.strip(".txt")

call_with_log(["mkdir", "-p fastqc"])
call_with_log(["mkdir", "-p " + projectname + "_vcf"])
call_with_log(["mkdir", "-p trim"])
call_with_log(["mkdir", "-p intermediate_files"])

filename = 'failed_'+InFileName
failedF = open(filename,'w')

filename2 = 'not_illumina_'+InFileName
failedplatform = open(filename2, 'w')

for line in InFile:
    numBP = 0
    line = line.strip()
    inputList = line.split()
    RGID = inputList[0]
    RGSM = inputList[1]
    RGLB = inputList[2]
    RGPL = inputList[3]
    pair = inputList[4]
    
    fastq_dump()
    print("fastq_dump completed")
    if RGPL != "illumina":
        failedplatform.write(line + '\n')
        failedF.write(RGID+ ' - not illumina' + '\n')
        continue
    if os.path.exists(RGID+'_2.fastq') != True and pair == 'paired':
        call_with_log(["rm", RGID+"_1.fastq"])
        failedF.write(RGID + ' - fastq_dump' + '\n')
    else:
        print('2 files were created')
        numBP = fastqc()
        print("fastqc completed")
        trim_galore()
        print('trim_galore completed')
        if numBP >= 70:
            bwaMEM()
        else:
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
            failedF.write(RGID + ' - bwa' + '\n')

InFile.close()
failedF.close()
failedplatform.close()

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

