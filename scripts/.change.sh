#!/bin/bash

# get counts
cat samples_adapters.txt | parallel -I ,_, "samtools view -F 4 ,_,.bam | cut -f 10,2 | sed 's:\t:_:g' | awk '{ cnts[\$0] += 1 } END { for (v in cnts) print cnts[v], v }' | sed -E 's/^ *([0-9]*) (.+)/\2\t\1/g' >> ,_,.tsv"

# find all unique sources
rm temp_builder.tsv

for f in *.bam; do
    samtools view -F 4 $f | cut -f 3,4,6 >> temp_builder.tsv
done

cat temp_builder.tsv | sort | uniq > temp_builder_2.tsv

rm temp_builder.tsv