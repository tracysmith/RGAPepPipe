#!/bin/bash

java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa -I ERR867518.realn.bam -o ERR867518.intervals

java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T IndelRealigner -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa -I ERR867518.realn.bam -targetIntervals ERR867518.intervals -nWayOut ERR867518.patrealn.bam

tar -zcvf ERR867518.realign.tar.gz realign/
