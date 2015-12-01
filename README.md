RGAPepPipe
==========

Reference Guided Assembly Pepperell Lab Pipeline

File Descriptions
==================
SRAtoVCF.py
------------
This script takes a list of SRA filenames (w/o .sra extension) and outputs a VCF file for each. For each sample, a tab-delimited .txt file is required containing:

> SRAReadGroupID	SampleGroupID	LibraryID	SeqPlatform	paired/single

For example:

> ERR07108	ERS088906	ERX048850	illumina	paired

The script must be run from a folder with downloaded .sra files. The path to the indexed reference .fa file must also be provided.

When running the script, type:
> python [path to SRAtoVCF.py] [path to .txt file] [.fasta reference file] [options]

For other available command line options run "python SRAtoVCF.py -h"

####Necessary programs required to run script:

fastq-dump
  * part of the SRAtoolkit software package available from NCBI
  * we use version 2.5.2
  * converts SRA format files to fastq  

fastqc
  * available from Babraham Bioinformatics, package includes necessary Picard BAM/SAM libraries
  * we use v0.11.3
  * performs quality checks on read data in fastq files  

trim-galore
  * available from Babraham Bioinformatics, fastqc and cutadapt are required
  * we use version 0.4.0 
  * trims bases below a given quality threshold, removes adapter sequences, and reruns fastqc on trimmed data  

bwa
  * available from bio-bwa.sourceforge.net
  * we use version 0.7.12 
  * maps reads to reference sequence  

SAMtools
  * available from samtools.sourceforge.net 
  * we use version 1.2 
  * sorts bam files by coordinates  

Picard-tools
  * available from picard.sourceforge.net
  * we use version 1.138
  * removes duplicates and edits read group information  

GATK (genome analysis toolkit)
  * available from broad institute
  * we use version 3.4.46
  * locally realigns reads and produces VCF  

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

VCFtoFasta.py
-------------
This script takes a variable number of input vcf files. Each .vcf file must contain 'confident sites' about one strain. At least one input file must be supplied. The name of the output file(s) will be automatically generated from the sample name in the vcf (SRS or ERS number) plus the identifying information contained in the name of the input file (ideally the RGID). This script is preliminary and has not been extensively tested. The script has been modified to produce a single contig fasta with gaps ("-"). For the original version that makes a multi-contig fasta look at # lines (new contig for every gap). Further devolpment should address a minimum size of contigs to be output (no longer relevant), how to deal with INDELS, and biopython fasta writing - currently entire sequence is output on one line.

> VCFtoFasta.py [input vcf file 1] .... [input vcf file n]

enaFileParser.py
----------------
This script converts the text file downloaded from the European Nucleotide Archive for a study and converts it into the input file for SRAtoVCF.py and a list of URLs of fastqs to be downloaded from ENA.

Usage:

> enaFileParser.py [input ENA text file]

get_url.py
---------
This script takes as input the ERPXXXXXX_download.txt file output by the enaFileParser.py and runs wget to download the fastq files. A '-t' flag allows the user to specify the number of threads.

Example
--------


In the template directory lies all the reference genome files for M. tuberculosis (the H37Rv files). Also included is a text file ERP000132.txt These are two of the three things necessary in order to run SRAtoVCF.py. The third file is the SRA file which can be obtained from NCBI using the following command in the command line:

> wget -m ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR027/ERR027082

This will create a folder of several sub directories in which ERR027082.sra is located. Move the .sra file to working directory. Then run the following command to execute SRAtoVCF.py:

> python SRAtoVCF.py ERP000132.txt H37Rv.fasta 

This will create a folder called ERP000132_vcf and within it will be a file ERR027082.vcf. You will use this file to run VCFtoSnpTable.py by running the command:

> python VCFsToSnpTable.py ERR027082.vcf SNPTable.txt

This will create a SNPTable.txt file which you will use to run SNPTableToNexus.py which is run with the following command:

> python SNPTableToNexus.py SNPTable.txt full.nex H37Rv.fasta

This will create full.nex which is the complete confident site nexus file. In order to get the SNP alignment nexus file run:

> python FulltoSNP.py full.nex end.nex .25

If all ran smoothly the end.nex file should match the final.nex file that is included in the template folder. 

