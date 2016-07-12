#!/bin/bash

java -ea -Xmx7g -jar /opt/PepPrograms/popoolation2_1201/mpileup2sync.jar --input ERR867518_Q20.mpileup --output ERR867518_Q20.sync --fastq-type sanger --min-qual 20 

tar -zcvf ERR867518.sync.tar.gz sync/
