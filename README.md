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

The script must be run from a folder with downloaded .sra files. There must also be a subdirectory called "trimfastqc" as well as a subdirectory titled "trim". The path to the indexed reference .fa file must also be provided.

When running the script, type:
> ipython [path to SRAtoVCF.ipy] [path to .txt file] [.fasta reference file] [path to directory named "trim"]

Necessary programs required to run script:

fastq-dump --> part of the SRAtoolkit software package available from NCBI (version we use 2.2.0).  
fastqc --> available from Babraham Bioinformatics, package includes necessary Picard BAM/SAM libraries (Version we use v0.10.1).  
trim-galore --> available from Babraham Bioinformatics, fastqc is required (we use version 0.2.6).  
bwa --> available from bio-bwa.sourceforge.net (we use version 0.7.7).  
SAMtools --> available from samtools.sourceforge.net (we use version 0.1.18).  
Picard-tools --> available from picard.sourceforge.net (we use version 1.107).  
GATK (genome analysis toolkit) --> available from broad institute (we use version 2.8.1).  

Java 7 must also be installed. 

VCFsToSnpTable.py
-----------------

This script takes a variable number of input vcf files. Each .vcf file must contain snp information about only one strain. At least one input file must be supplied. The name of an output file must also be supplied, and output is in snp table format. The output contains the combined information about all strains. 

> VCFsToSnpTable.py [input .vcf 1] .... [input.vcf n] [output file]


SNPTableToNexus.py
-----------------

Takes the SnpTable generated in VCFsToSnpTable.py and converts it to a nexus file. Requires a reference genome as well. Output file should end in .nex.   

Usage:  

> SNPTableToNexus.py [input file] [outputfile] [reference sequence .fa file]

FulltoSNP.py
------------
Takes the nexus file generated in SNPTableToNexus.py and performs a SNP alignment. Output file should be in .nex format. User must also input a threshold value (between 0 and 1) . The threshold determines the number of sequences with an ambiguous base at a certain position and rejects base position if the percentage of ambiguous bases is greater than the threshold.  
USAGE:

> FulltoSNP.py [input nexus file] [output file] [threshold value] 


Example
--------


In the template directory lies all the reference genome files for M.Tb (the H37Rv files). Also included is a text file ERP000132.txt These are two of the three things necessary in order to run SRAtoVCF.ipy. The third file is the SRA file which can be obtained from NCBI using the following command in the command line:

> wget -m ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR027/ERR027082

This will create a folder of several sub directories in which ERR027082.sra is located. Move the .sra file to working directory. There should also be to folders in the current working directory that you need to make. trim and trimfastqc , Then run the following command to execute SRAtoVCF.ipy:

> ipython SRAtoVCF.ipy ERP000132.txt H37Rv.fasta ./trim

This will create a folder called ERP000132_vcf and within it will be a file ERR027082.vcf. You will use this file to run VCFtoSnpTable.py by running the command:

> python VCFsToSnpTable.py ERR027082.vcf SNPTable.txt

This will create a SNPTable.txt file which you will use to run SNPTableToNexus.py which is run with the following command:

> python SNPTableToNexus.py SNPTable.txt full.nex H37Rv.fasta

This will create full.nex which is the complete confident site nexus file. In order to get the SNP alignment nexus file run:

> python FulltoSNP.py full.nex end.nex .25

If all ran smoothly the end.nex file should match the final.nex file that is included in the template folder. 

