#!/bin/bash

java -Xmx2g -jar /opt/PepPrograms/GenomeAnalysisTK.jar -I $1.ready.bam -R $2 -T RealignerTargetCreator -o $1.intervals

java -Xmx2g -jar /opt/PepPrograms/GenomeAnalysisTK.jar -I $1.ready.bam -R $2 -T IndelRealigner -targetIntervals $1.intervals -o $1.realn.bam
