#!/bin/bash
bed=$1
genome=$2
name=$3
#summit expansion to 51bp
bedtools slop -i ${bed} -g ${genome} -b 25 | awk -F "\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' > ${name}_summits_50bp.tab
#summit expansion to 401bp
bedtools slop -i ${bed} -g ${genome} -b 200 | awk -F "\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' > ${name}_summits_400bp.tab
