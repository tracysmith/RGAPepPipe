#!/bin/bash

file1=`find $1 -type f -print`

if [ -n "$file1" ]; then #n flag means if string is NOT null -z means IS null
    exit 1
fi
