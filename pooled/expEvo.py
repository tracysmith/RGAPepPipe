#!/usr/bin/python

from subprocess import call

#3.9.2015 - Writing a script to process Tracy's experimentally evolved
#Mtb samples.

#There are 6 populations that were evolved, with a total of 18 samples.
samples = [
    "ERR003100", "ERR003108", "ERR003112", "ERR003116",
    "ERR004900", "ERR004908", "ERR004912", "ERR004916",
    "ERR005500", "ERR005504",
    "ERR007200", "ERR007204", "ERR007208",
    "ERR034500", "ERR034504", "ERR034508", 
    "ERR054000", "ERR054004"
    ]
print("There are " + str(len(samples)) + " samples in the list.") 

def remove_ambMap(samp):
    #Remove ambiguously mapped reads from the realigned bam files 
    call('samtools view -q 20 -b /home/tmsmith/data/expEvo/RGA/RGAbamsBais/{samp}.realn.bam' '| samtools sort - {samp}'.format(samp=samp), shell=True)

#Convert the relevant bam files into an mpileup (by population, 
#i.e. all samples from #31 together):
#samtools mpileup -B -f ./<path to ref>/MtbNCBIH37Rv.fa -q 20 \
#-Q 20 ./<path to bams - sep by space>/ > {Prefix}.pileup 

#Convert the mpileups to sync files:
#java -ea -Xmx7g -jar \
#/home/peplab/src/popoolation2_12-1/mpileup2syn.jar \
#--input {Prefix}.mpileup --output {Prefix}.sync --fastq-type sanger \
#--min-qual 20 --threads 4

def callable_loci(samp):
    """Run callable loci on samples"""
   # call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I /home/tmsmith/data/expEvo/RGA/RGAbamsBais/{samp}.realn.bam -summary {samp}_defaults.summary -o {samp}_defaults.bed -R /home/mrood/data/Mtb/ref/MtbNCBIH37Rv.fa'.format(samp=samp), shell=True)
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I /home/tmsmith/data/expEvo/RGA/RGAbamsBais/{samp}.realn.bam -summary {samp}_strict.summary -o {samp}_strict.bed -R /home/mrood/data/Mtb/ref/MtbNCBIH37Rv.fa -frlmq 0.04 -mmq 20'.format(samp=samp), shell=True)


#Prefix = raw_input("Prefix: ")
def snp_frequency_diff(Prefix):
    #Calculate allele frequency differences
    call('perl /opt/PepPrograms/popoolation2_1201/snp-frequency-diff.pl --input {Prefix}.sync --output {Prefix}_mc10 --min-count 10 --min-coverage 10 --max-coverage 2%'.format(Prefix=Prefix), shell=True)

def fisher_test(Prefix):
    #Estimate the significance of allele frequency differences
    call('perl /opt/PepPrograms/popoolation2_1201/fisher-test.pl --input {Prefix}.sync --output {Prefix}_mc10.fet --min-count 10 --min-coverage 10 --max-coverage 2% --min-covered-fraction 1 --window-size 1 --step-size 1 --suppress-noninformative'.format(Prefix=Prefix), shell=True)

def fst_sliding(Prefix):
    #Calculate Fst values using a sliding-window approach
    call('perl /opt/PepPrograms/popoolation2_1201/fst-sliding.pl --input {Prefix}.sync --output {Prefix}_mc10_p10K.fst --min-count 10 --min-coverage 10 --max-coverage 2% --min-covered-fraction 1 --window-size 1 --step-size 1 --suppress-noninformative --pool-size 10000'.format(Prefix=Prefix), shell=True)

for samp in samples:
    print("Processing sample: " + samp)
#    remove_ambMap(samp)
    callable_loci(samp)

#snp_frequency_diff(Prefix)
#fisher_test(Prefix)
#fst_sliding(Prefix)

def make_fetDict(Prefix):
    fetDict = {}
    with open('{Prefix}_mc10.fet'.format(Prefix=Prefix), 'r') as fetFile:
        for line in fetFile:
            line = line.strip().split('\t')
            pos = int(line[1])
            fetDict[pos] = {}
            paircomps = line[5:]
            for i in paircomps:
                key,fet = i.split("=")
                fetDict[pos][str(key)] = float(fet) 
    return fetDict

def make_figureFile(Prefix, fetDict):
    with open("{Prefix}_mc10_p10K.fst".format(Prefix=Prefix), 'r') as fstFile, open("{Prefix}_mc10_manhatten.txt".format(Prefix=Prefix), 'w') as outFile:
        outFile.write("%s\t%s\t%s\t%s\n" %
        ("key",
        "position",
        "fst",
        "fet")
        )
        for line in fstFile:
            line = line.strip().split('\t')
            pos = int(line[1])
            paircomps = line[5:]
            for i in paircomps:
                key = str(i.split("=")[0])
                fst = float(i.split("=")[1])
                outFile.write("%s\t%i\t%f\t%f\n" %
                (key,
                pos,
                fst,
                fetDict[pos][key])
                )

#fetDict = make_fetDict(Prefix)
#make_figureFile(Prefix, fetDict)
