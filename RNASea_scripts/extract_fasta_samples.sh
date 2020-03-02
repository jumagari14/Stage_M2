#!/bin/bash

lines=$(($1*$2))


currentD=$(readlink -f $3)

mkdir -p -m 755 "$currentD"/test_data

list=$(find $currentD -name "*.fastq.*")

for i in $list 
do 
    filename=$(echo "$i" | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1)
    gzip -dc "$i" | cat | head -$lines > "$currentD"/test_data/"$filename".fq
done 

exit 0 