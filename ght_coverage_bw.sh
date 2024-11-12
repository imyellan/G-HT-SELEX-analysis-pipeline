#!/bin/bash

outdir=~/TEHMM_proj/analysis/selex_results/ght_explore
mkdir -p "$outdir"
in_tab=$1
gene_nm="$(basename "$(dirname "$in_tab")")"
gene_out="$outdir/$gene_nm"

mkdir -p "${gene_out}"
Rscript ~/TEHMM_proj/selex_motif_scripts/make_coverage_bedgraphs.R "$in_tab" \
"${gene_out}"

for bg in "${gene_out}"/*.bed; do
    bedGraphToBigWig $bg ~/data/hg38.chrom.sizes "${bg%%.bed}.bw"
done
rm "${gene_out}"/*.bed