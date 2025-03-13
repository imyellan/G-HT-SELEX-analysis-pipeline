#!/usr/bin/env bash
dtime=$(date +"%y%m%d%H%M")
# cd $HOME
export NXF_SINGULARITY_HOME_MOUNT=true ## required for singularity to work properly

submitjob -m 1 -w 20 -E isaac.yellan95@gmail.com \
cd "$HOME"\; nextflow  -c ~/TEHMM_proj/selex_motif_scripts/nextflow.config \
run ~/TEHMM_proj/selex_motif_scripts/motif_discovery_eval.nf \
--inputDir "/home/hughespub/SELEX_Data/GHT-SELEX_Fastqs_Ready_For_SRA_etc/YWWX_250205_A00546_0203_AHT2JVDRX5_HUGHES" \
--inputType "fastq" \
--inputExp "GHT" \
--outputDir "${HOME}/TEHMM_proj/analysis/selex_results/motif_pipeline" \
-resume \
-with-report ~/nextflow_reports/"${dtime}"_fl.html -N isaac.yellan95@gmail.com \
-with-dag ~/nextflow_reports/"${dtime}"_dag.html \
-with-timeline ~/nextflow_reports/"${dtime}"_timeline.html \
-qs 999 \
-profile cluster -ansi-log false \&\> ~/nextflow_log.txt

# nextflow run ~/TEHMM_proj/selex_motif_scripts/motif_discovery_eval.nf \
# --inputDir "/home/hughespub/SELEX_Data/GHT-SELEX_Fastqs_Ready_For_SRA_etc/YWWX_250205_A00546_0203_AHT2JVDRX5_HUGHES" \
# --inputType "fastq" \
# --inputExp "GHT" \
# --outputDir "${HOME}/TEHMM_proj/analysis/selex_results/motif_pipeline"-resume \
# -profile local


# params.ywuv_htselex = "/home/hughespub/SELEX_Data/HT-SELEX_Fastqs_Ready_For_SRA_etc/YWUV_and_RoziHadiAttackATAC_240711_A00546_0178_AH57KMDRX5"
# params.ywt_htselex = "/home/hughespub/finishTF/YWT_Mostly_Isaac/FinalFastq_Trimmed_and_Filtered_Ready"
# params.yww_htselex = "/home/hughespub/SELEX_Data/HT-SELEX_Fastqs_Ready_For_SRA_etc/HTSELEX_BcorrectFilteredTrimmed"
# // params.ght_input1 = "/home/hughespub/SELEX_Data/GHT-SELEX_Fastqs_Ready_For_SRA_etc/YWUV_and_RoziHadiAttackATAC_240711_A00546_0178_AH57KMDRX5"
# params.ght_input = "/home/hughespub/SELEX_Data/GHT-SELEX_Fastqs_Ready_For_SRA_etc"
# params.ght_peaks = "/home/hughespub/ahcorcha/Transposone_TFs_MAGIX_peaks_08_11_24"