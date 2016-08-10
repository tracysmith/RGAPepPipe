#!/usr/bin/env python

import os
import sys

# take input file and check for pass fail 
fastqc_1 = open(sys.argv[1],'r')
fastqc_2 = open(sys.argv[2],'r')

run = sys.argv[1].split("_")[0]

# initializign results as False, will make true once pass
result_1 = False
result_2 = False

def checkQual(fastqc):
    results = []
    for line in fastqc:
        line = line.strip().split("\t")
        results.append(line[0])
    print results
    numFails = 0
    for result in results:
        if result == "FAIL":
            numFails += 1
    print numFails
    if numFails >= 5:
        return False
    else:
        return True

result_1 = checkQual(fastqc_1)
result_2 = checkQual(fastqc_2)

if result_1 and result_2 == True:
    print "whoo all clear"
else: # write failed file
    outfileName = str(run) + "_failedFastqc.txt"
    outfile = open(outfileName,'w')
    outfile.write("The fastqc results are below the specified threshold.\n")
    print "the dag is done-zo"

