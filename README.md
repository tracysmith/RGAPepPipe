RGAPepPipe
==========

Reference Guided Assembly Pepperell Lab Pipeline

File Descriptions
==================
SRAtoVCF.ipy
------------
This script takes a list of SRA filenames (w/o .sra extension) and outputs a VCFfile for each. For each sample, a tab-delimited .txt file is required containing:

> SRAReadGroupID	SampleGroupID	LibraryID	SeqPlatform	paired/single

For example:

> ERR07108	ERS088906	ERX048850	illumina	paired

The script must be run from folder with indexed reference genome files in it as well as downloaded .sra files. There must also be a subdirectory called "trimfastqc" as well as a subdirectory titled "trim"

When running the script, type:
> ipython [path to SRAtoVCF.ipy] [path to .txt file] [.fasta reference file] [path to directory named "trim"]

Necessary programs required to run script:

fastq-dump --> part of the SRAtoolkit software package available from NCBI (version we use 2.2.0).  
fastqc --> available from Babraham Bioinformatics, package includes necessary Picard BAM/SAM libraries (Version we use v0.10.1).  
trim-galore --> available from Babraham Bioinformatics, fastqc is required (we use version 0.2.6).  
bwa --> available from bio-bwa.sourceforge.net (we use version 0.6.1-r104).  
SAMtools --> available from samtools.sourceforge.net (we use version 0.1.18).
Picard-tools --> available from picard.sourceforge.net (we use version 1.107).  
GATK (genome analysis toolkit) --> available from broad institute (we use version 2.8.1).  

Java must also be installed. 

VCFsToSnpTable.py
-----------------

blah


