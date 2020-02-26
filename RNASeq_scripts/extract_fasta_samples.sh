#!/bin/bash

currentD=$(pwd)

mkdir -p -m 755 Stage\ M2/test_data

list=$(find $currentD -name "*.fastq.*")

for i in $list 
do 
    filename=$(echo "$i" | rev | cut -d "/" -f 1 | rev)
    zcat "$i" | head -$1 > Stage\ M2/test_data/"$filename"
done 

exit 0 