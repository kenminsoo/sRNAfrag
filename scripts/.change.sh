#!/bin/bash

cat samples_adapters.txt | parallel -I ,_, "samtools view -F 4 ,_,.bam | cut -f 10,2 | sed 's:\t:_:g' | awk '{ cnts[\$0] += 1 } END { for (v in cnts) print cnts[v], v }' | sed -E 's/^ *([0-9]*) (.+)/\2\t\1/g' >> ,_,.tsv"