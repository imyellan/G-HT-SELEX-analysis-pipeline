#!/usr/bin/env bash
dtime=$(date +"%y%m%d%H%M")
# cd $HOME
export NXF_SINGULARITY_HOME_MOUNT=true ## required for singularity to work properly

submitjob -m 1 -w 72 -E isaac.yellan95@gmail.com \
cd "$HOME"\; nextflow run ~/TEHMM_proj/selex_motif_scripts/motif_discovery_eval.nf -resume \
-with-report ~/nextflow_reports/"${dtime}"_fl.html -N isaac.yellan95@gmail.com \
-profile cluster -ansi-log false \&\> ~/nextflow_log.txt

# nextflow run ~/TEHMM_proj/selex_motif_scripts/motif_discovery_eval.nf -resume -profile local
