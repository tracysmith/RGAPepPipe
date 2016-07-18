#!/bin/bash

java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $1 "${@:3}" -o $2.intervals

java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T IndelRealigner -R $1 "${@:3}" -targetIntervals $2.intervals -nWayOut .patrealn.bam
