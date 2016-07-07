#!/bin/bash

perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input ERR867518_Q20.pileup --output ERR867518_Q20.indelreg.gtf

perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input ERR867518_Q20.pileup --gtf ERR867518_Q20.indelreg.gtf --output ERR867518_Q20.noIndel.pileup

tar -zcvf ERR867518.remIndels.tar.gz remIndels/
