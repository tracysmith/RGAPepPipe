#!/bin/bash

java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I ERR867518.realn.patrealn.bam -summary ERR867518_defaults.summary -o ERR867518_defaults.bed -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa 

java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I ERR867518.realn.patrealn.bam -summary ERR867518_strict.summary -o ERR867518_strict.bed -R /opt/data/mtuberculosis/MtbNCBIH37Rv.fa -frlmq 0.04 -mmq 20

tar -zcvf ERR867518.callableLoci.tar.gz callableLoci/
