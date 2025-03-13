#!/usr/bin/env bash

TF=$1
motif_pfm=$2
peak_bed=$3
outdir=~/TEHMM_proj/analysis/selex_results/ht_ght_comparison
mkdir -p "$outdir"

# echo $TF $motif_pfm $peak_bed
# peak_bed=~/TEHMM_proj/analysis/selex_results/motif_pipeline/target_C17ORF113_LTR_results_all_with_eFDR_PEAKS_kneedle/knee_filt_peaks.bed
# motif_pfm=~/TEHMM_proj/analysis/selex_results/motif_pipeline/YWU_A_3_C08_pTH14640_CC40NCCAATA_eGFP_IVT_S434_R1_001/meme_out_8/meme_pfm_1.pfm

# convert motif pfm to format compatible with moods
~/scripts/pwm_transpose_simple.sh "${motif_pfm}" > "${motif_pfm%.pfm}_moods.pfm"

# scan hg38 with motif, using p-value threshold of 0.001 (same as Codebook)
moods-dna.py -m "${motif_pfm%.pfm}_moods.pfm" -s ~/data/hg38_main_chr.fa -p 0.000001 \
--log-base 2 > "${motif_pfm%.pfm}_moods_hg38_000001.csv"

# convert moods output to bed format
awk -F"," 'BEGIN{OFS="\t"}{print $1,$3,$3+length($6),$2,$5,$4}' \
"${motif_pfm%.pfm}_moods_hg38_000001.csv" > "${motif_pfm%.pfm}_moods_hg38_000001.bed"
motif_bed="${motif_pfm%.pfm}_moods_hg38_000001.bed"

# intersect motif hits with peak bed
bedtools intersect -a "${peak_bed}" -b "${motif_bed}" -wao -F 1 \
> "${outdir}"/"${TF}_ht_motif_peak_intersect.txt"
# also run bedtools fisher test
bedtools fisher -a <(sort -k 1,1 -k2,2n "${peak_bed}") \
-b <(sort -k 1,1 -k2,2n "${motif_bed}") \
-g ~/data/hg38.chrom.sizes.sorted -F 1 \
> "${outdir}"/"${TF}_ht_motif_peak_fisher.txt"

## monte carlo simulation method to calculate enrichment fold change & emprical p-value
observed=$(bedtools intersect -a "${peak_bed}" -b "${motif_bed}" -F 1 -u | wc -l)
echo "Observed Overlap: $observed" > "${outdir}"/"${TF}_ht_motif_peak_enrichment_stats.txt"

run_shuffle () {
    peak_bed=$1
    motif_bed=$2

    bedtools shuffle -i "${peak_bed}" -g ~/data/hg38.chrom.sizes.sorted \
    | bedtools intersect -a - -b "${motif_bed}" -F 1 -u | wc -l
}
export -f run_shuffle

iterations=1000
total=0
shuffle_counts=()

for i in $(seq 1 $iterations); do
    shuffle_counts+=($(run_shuffle "${peak_bed}" "${motif_bed}"))
done

mapfile -t shuffle_counts < <(parallel -j 4 run_shuffle "${peak_bed}" "${motif_bed}" ::: $(seq 1 $iterations))
total=$(echo "${shuffle_counts[@]}" | tr ' ' '\n' | awk '{s+=$1} END{print s}')
expected=$(echo "$total / $iterations" | bc -l)

# calculate fold change
echo "Expected Overlap (Mean of $iterations shuffles): $expected" \
>> "${outdir}"/"${TF}_ht_motif_peak_enrichment_stats.txt"
foldchange=$(echo "$observed / $expected" | bc -l)
echo "Fold Change: $foldchange" >> "${outdir}"/"${TF}_ht_motif_peak_enrichment_stats.txt"

# calculate empirical p-value
greater_equal=$(echo "${shuffle_counts[@]}" | tr ' ' '\n' | awk -v obs=$observed '$1 >= obs' | wc -l)
p_value=$(echo "$greater_equal / $iterations" | bc -l)
echo "Empirical p-value: $p_value" >> "${outdir}"/"${TF}_ht_motif_peak_enrichment_stats.txt"
