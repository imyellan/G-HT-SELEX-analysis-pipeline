#!/usr/bin/env bash
micromamba deactivate
# cd $HOME
submitjob -m 1 -w 1 -E isaac.yellan95@gmail.com \
nextflow run ~/TEHMM_proj/selex_motif_scripts/motif_discovery_eval.nf -resume \
\&\> ~/nextflow_log.txt