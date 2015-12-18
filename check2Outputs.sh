#!/bin/bash

file1=`find $1 -type f -size +$3 -print`
file2=`find $2 -type f -size +$3 -print`

if [ -z "$file1" ]; then
    exit 1
fi

if [ -z "$file2" ]; then
    exit 1
fi
