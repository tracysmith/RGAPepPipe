#!/bin/bash

java -jar /opt/PepPrograms/pilon-1.16.jar --genome $2 --frags $1.realn.bam --output $1_pilon --variant --mindepth 10 --minmq 40 --minqual 20
