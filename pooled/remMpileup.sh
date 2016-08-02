#!/bin/bash

cat "${@:2}" > $1.indelreg.gtf

perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input $1_Q20_q20.mpileup --gtf $1.indelreg.gtf --output $1.noIndel.mpileup
