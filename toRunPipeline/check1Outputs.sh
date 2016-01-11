#!/bin/bash

file1=`find $1 -type f -size +$2 -print`

if [ -z "$file1" ]; then
    exit 1
fi
