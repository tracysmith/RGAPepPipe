#!/bin/bash

perl /opt/PepPrograms/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ERR867518_Q20.mpileup --output ERR867518.indelreg.gtf

perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input ERR867518_Q20.mpileup --gtf ERR867518.indelreg.gtf --output ERR867518.noIndel.mpileup

tar -zcvf ERR867518.remMpileup.tar.gz remMpileup/
