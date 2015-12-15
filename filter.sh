#!/bin/bash

java -Xmx2g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T VariantFiltration -R $2 -o $1_filter.vcf --variant $1.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "RGAPepPipeFilter"
