#!/usr/bin/env bash

cycle=$1
fastq=$2

reformat.sh in=${fastq} out=cycle_${cycle}.fa
$HOME/scripts/fasta_linearize.sh cycle_${cycle}.fa \
| shuf -n 190000 \
| awk -v c=$cycle 'BEGIN{OFS="\n"} {print "> peak: 20; signal: "c, $2}'
