#!/bin/env bash

proj_dir=/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline
n_jobs=20
rm /home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/selex_res_motif_bmark_*
num_bmarks=$(find "${proj_dir}" -name "*eval_summ.txt" | wc -l)
recalc=FALSE

for i in "${proj_dir}"/*/pos.fa.gz; do
    zgrep ">" "$i" | wc -l | sed -E 's/ +//g' > "${i%%.fa.gz}"_seq_count.txt
done

# create benchmark result summary df
Rscript /home/hugheslab1/iyellan/TEHMM_proj/selex_motif_scripts/generate_bmark_summ_df.R TRUE

## divide the number of benchmarks by the number of jobs, create a list of ranges
## to be used for parallel processing
start_idx=(1)
n_bmarks_per_job=$((num_bmarks/n_jobs))
end_idx=($((start_idx[0] + n_bmarks_per_job)))

# create ranges for parallel processing
for i in $(seq 0 $((n_jobs-1))); do
    if [[ "$i" -gt 0 ]]; then
        start_idx[i]=$((end_idx[$((i-1))]+1))
        end_idx[i]=$((start_idx[i]+n_bmarks_per_job))
    fi
    if [[ ${end_idx[i]} -gt $num_bmarks ]]; then
        end_idx[i]=$num_bmarks
    fi
    if [[ "$recalc" == TRUE ]]; then
        mem=75
        threads=5
        wtime=5
    else
        mem=5
        threads=5
        wtime=1
    fi
    submitjob -w $wtime -m $mem -c $threads -E isaac.yellan95@gmail.com \
    Rscript /home/hugheslab1/iyellan/TEHMM_proj/selex_motif_scripts/motif_analysis_parse.R \
    ${start_idx[i]} ${end_idx[i]} $recalc
done

