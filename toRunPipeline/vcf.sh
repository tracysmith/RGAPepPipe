#!/bin/bash

java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -I $1.realn.bam -R $2 -T UnifiedGenotyper -o $1.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -stand_call_conf 20 -stand_emit_conf 20 --sample_ploidy 1 -nt 4 -rf BadCigar
