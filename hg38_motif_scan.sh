#!/bin/env bash

in_pfm=$1
peak_bed=$2

## transpose for compatibility with moods
awk '
BEGIN{FS=OFS="\t"}
{
    for (i = 1; i <= NF; i++)  {
        a[i, NR] = $i
    }
}
{
    for (i = 1; i <= NF; i++)  {
        a[i, NR] = $i
    }
}
NF > p { p = NF }
END {
    for (i = 1; i <= p; i++) {
        for (j = 1; j <= NR; j++) {
            printf "%s%s", a[i, j], (j == NR ? "\n" : " ")
        }
    }
}' "$in_pfm" > "${in_pfm%.pfm}"_transposed.pfm

moods-dna.py -m "${in_pfm%.pfm}"_transposed.pfm -s /home/hugheslab1/iyellan/C17ORF113_all_peaks.fa\
-p 0.00001 -o "${in_pfm%.pfm}"_hg38_moods_out.txt --ps 0.0001 --log-base 10

rm "${in_pfm%.pfm}"_transposed.pfm

# output bed file version of the moods output
awk -F"," 'BEGIN{OFS="\t"} {print $1,$3,length($6)+$3,$2,$5,$4}' "${in_pfm%.pfm}"_hg38_moods_out.txt \
> "${in_pfm%.pfm}"_hg38_moods_out.bed

# bedtools intersect to get the hg38 coordinates
bedtools intersect \
-a <(tail -n +2 ${peak_bed}) \
-b "${in_pfm%.pfm}"_hg38_moods_out.bed -F 1 -wao \
> "${in_pfm%.pfm}"_moods_MAGIX_intersect.tsv

gzip -f "${in_pfm%.pfm}"_moods_MAGIX_intersect.tsv "${in_pfm%.pfm}"_hg38_moods_out.bed \
"${in_pfm%.pfm}"_hg38_moods_out.txt