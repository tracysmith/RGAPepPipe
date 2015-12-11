#!/bin/bash

java -Xmx2g -jar /opt/PepPrograms/picard-tools-1.138/picard.jar MarkDuplicates I=$1.sort.bam O=$1.dedup.bam M=$1.metrics REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT
java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/picard.jar AddOrReplaceReadGroups I=$1.dedup.bam O=$1.ready.bam RGID=$1 RGLB=$2 RGPL=$3 RGPU=dummy-barcode RGSM=$4 VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true
