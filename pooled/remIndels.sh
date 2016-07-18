#!/bin/bash

perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input $1.mpileup --output $1.indelreg.gtf

perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $1.mpileup --gtf $1.indelreg.gtf --output $1.noIndel.pileup
