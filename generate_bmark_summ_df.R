library(dplyr)
library(readr)
library(ggplot2)
library(purrr)

proj_dir <<- ifelse(Sys.info()["sysname"] == "Darwin",
                    "~/iyellan_bc2-clust/TEHMM_proj/analysis/selex_results/motif_pipeline",
                    "/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline")
reimport_bmarks <- ifelse(Sys.info()["sysname"] == "Darwin", 
                          F, as.logical(commandArgs(trailingOnly = T)[1]))

import_benchmark_summ <- function(bmark_path){
  pTH <- ifelse(grepl("PEAKS", bmark_path), 
                basename(dirname(bmark_path)),
                gsub(".*_(pTH[0-9]+)_.*", "\\1", bmark_path))
  motif <- gsub("_[0][0-9]+", "", gsub("_eval_summ.txt", "", basename(bmark_path)))
  bmark_pos_frac <- gsub("^0", "0.", 
                         gsub(".*_([0][0-9]+)_.*", "\\1", basename(bmark_path)))
  cycle <- ifelse(grepl("merged", bmark_path), "merged",
                  ifelse(grepl("PEAKS", bmark_path), "MAGIX", 
                         ifelse(grepl("Cycle", bmark_path), 
                                gsub(".*_Cycle([1-3]).*", "\\1", bmark_path),
                                gsub(".*_A_([1-3])_.*", "\\1", bmark_path))))
  well <- ifelse(grepl("PEAKS", bmark_path), "MAGIX",
                 gsub(".*_([A-H][0-1][0-9])_.*", "\\1", dirname(bmark_path)))
  selex_plate <- ifelse(grepl("PEAKS", bmark_path), "MAGIX",
                        gsub(".*(YW[TUV])_.*", "\\1", bmark_path))
  exp_id <- basename(dirname(bmark_path))
  
  bmark_df <- read_table(bmark_path, col_names = c("metric", "score"), progress = F,
                         col_types = c("c", "d")) %>%
    mutate(pTH = pTH, motif = motif, cycle = cycle, Well = well, selex_plate = selex_plate,
           bmark_pos_frac = bmark_pos_frac, exp_id = exp_id)
  bmark_df
}

## import AUROC and AUPRC values
if (reimport_bmarks) {
  benchmark_summ_fils <- list.files(
    proj_dir,
    pattern = "eval_summ.txt", full.names = T, recursive = T)
  benchmark_summ_df <- benchmark_summ_fils %>% map(import_benchmark_summ) %>%
    list_rbind()
  write_csv(benchmark_summ_df, file.path(proj_dir, "../motif_benchmark_summary_df.csv.gz"))
} else{
  benchmark_summ_df <- read_csv(file.path(proj_dir, "../motif_benchmark_summary_df.csv.gz"))
}


