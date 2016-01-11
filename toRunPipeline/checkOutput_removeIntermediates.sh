#!/bin/bash

file=`find $1 -type f -size +$2 -print`
if [ -z "$file" ]; then
    exit 1
fi

rm "${@:3}"
