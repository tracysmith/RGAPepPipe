#!/bin/bash

java -ea -Xmx7g -jar /opt/PepPrograms/popoolation2_1201/mpileup2sync.jar --input $1.noIndel.mpileup --output $1.sync --fastq-type sanger --min-qual 20 
