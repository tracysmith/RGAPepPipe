#!/bin/bash

unzip $1_1_val_1_fastqc.zip
unzip $1_2_val_2_fastqc.zip

python fastqc_check.py $1_1_val_1_fastqc/summary.txt $1_2_val_2_fastqc/summary.txt


