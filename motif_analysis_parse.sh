#!/bin/bash

proj_dir=/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline
n_jobs=25

## generate benchmark df
Rscript ~/TEHMM_proj/selex_motif_scripts/generate_bmark_summ_df.R TRUE \
2> /dev/null

num_bmarks=$(find "${proj_dir}" -name "*eval_summ.txt" | wc -l)
n_bmarks_per_job=$((num_bmarks/n_jobs))
## divide the number of benchmarks by the number of jobs, create a list of ranges
## to be used for parallel processing
start_arr_last=1
unset start_arr end_arr
declare -a start_arr
declare -a end_arr
# start_arr=(${start_arr_last})
# end_arr=($((sstart_arr_last + n_bmarks_per_job)))

while [[ $((start_arr_last + n_bmarks_per_job)) -lt ${num_bmarks} ]]; do
    start_arr+=(${start_arr_last})
    end_arr+=($((start_arr_last+n_bmarks_per_job)))
    start_arr_last=$((${end_arr[-1]}+1))
done
## add the last range
start_arr+=($((${end_arr[-1]}+1)))
end_arr+=($num_bmarks)


## run the parallel jobs
for i in $(seq 0 $((n_jobs-1))); do
    echo \
    Rscript ~/TEHMM_proj/selex_motif_scripts/motif_analysis_parse.R \
    "${start_arr[$i]}" "${end_arr[$i]}"
done > ~/motif_analysis_parse_jobs.sh

## break up the jobs into smaller chunks no larger than 10 lines
split -l 10 ~/motif_analysis_parse_jobs.sh ~/motif_analysis_parse_jobs_
rm ~/motif_analysis_parse_jobs.sh
for i in ~/motif_analysis_parse_jobs_*; do
    submitjob -f "$i" -m 2 -w 2 -E isaac.yellan95@gmail.com
done
