#!/bin/bash

perl /opt/PepPrograms/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input $1_Q20_q20.mpileup --output $1.indelreg.gtf

perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $1_Q20_q20.mpileup --gtf $1.indelreg.gtf --output $1.noIndel.mpileup
