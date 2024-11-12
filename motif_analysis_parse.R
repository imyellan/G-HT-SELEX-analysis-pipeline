library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(furrr)
library(pracma)
library(magrittr)

proj_dir <<- ifelse(Sys.info()["sysname"] == "Darwin",
                   "~/iyellan_bc2-clust/TEHMM_proj/analysis/selex_results/motif_pipeline",
                   "/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline")
range_start <- as.numeric(commandArgs(trailingOnly = T)[1])
range_end <- as.numeric(commandArgs(trailingOnly = T)[2])
benchmark_summ_df <- read_csv(file.path(proj_dir, "../motif_benchmark_summary_df.csv.gz"))

import_ROC <- function(f){
  read_tsv(f) %>% 
    mutate(base = gsub("_ROC.tsv", "", gsub(".*motif_pipeline/", "", f)))
}

min_max_normalization <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

## alternative normalization using custom boundaries; i.e. how much of a set
## area the curve integral takes up
custom_normalization <- function(x, cust_max) {
  return((x - 0) / (cust_max - 0))
}

spline_merge_pAUROC <- function(ROC_res_df, ROC_spline_df, thresh){
  ROC_res_df_thresh <- ROC_res_df %>% filter(fpr <= thresh)
  if (max(ROC_res_df_thresh$fpr) < thresh) {
    if(nrow(ROC_res_df_thresh) == 0){
      ROC_res_df_thresh <- data.frame(fpr = 0, tpr = 0)
    }
    ROC_spline_df %<>% filter(fpr <= thresh)
    if(nrow(ROC_spline_df) == 0){
      ROC_spline_df <- data.frame(fpr = 0, tpr = 0)
    }
    unique(rbind(ROC_res_df_thresh %>% select(fpr, tpr), 
                 ROC_spline_df %>% select(fpr, tpr))) %$%
      trapz(fpr, tpr)
  }
  else{
    trapz(ROC_res_df_thresh$fpr, ROC_res_df_thresh$tpr)
  }
}

part_auroc <- function(roc_path, recalc = F){
  pauroc_df_path <- gsub("_ROC.tsv", "_pauroc_df.csv.gz", file.path(proj_dir, roc_path))
  if (recalc | !file.exists(pauroc_df_path)) {
    ROC_res_df <- import_ROC(file.path(proj_dir, roc_path))
    ## while at it, save an ROC plot
    ROC_res_df %>% 
      ggplot(aes(fpr, tpr)) + geom_line() +
      ggtitle(roc_path) + 
      theme_bw() + theme(legend.position = "none")
    ggsave(gsub(".tsv", ".pdf", file.path(proj_dir, roc_path)))
    
    ## create a spline of the ROC curve
    ROC_spline <- spline(x = ROC_res_df$fpr, ROC_res_df$tpr, 
                         n = 1000*length(ROC_res_df$fpr), ties = "max")
    ROC_spline_df <- data.frame(fpr = ROC_spline$x, tpr = ROC_spline$y)
    
    ## save partial ROC plot - using FPR <= 1e-3, TPR <= 1e-2
    ROC_res_df %>% 
      ggplot(aes(fpr, tpr)) + geom_line() + 
      ggtitle(roc_path) + 
      theme_minimal() + theme(legend.position = "none") +
      coord_cartesian(xlim = c(0, 1e-3), ylim = c(0, 1e-2))
    ggsave(gsub(".tsv", "_partial.pdf", file.path(proj_dir, roc_path)))
    
    # calculate normalized pAUROC at FPR <= 1e-3
    # part_auroc_1e3 <- trapz(roc_tbl_fpr_1e3$fpr, roc_tbl_fpr_1e3$tpr)
    # part_auroc_1e4 <- trapz(roc_tbl_fpr_1e4$fpr, roc_tbl_fpr_1e4$tpr)
    # part_auroc_1e5 <- trapz(roc_tbl_fpr_1e5$fpr, roc_tbl_fpr_1e5$tpr)
    # part_auroc_1e10 <- trapz(roc_tbl_fpr_1e10$fpr, roc_tbl_fpr_1e10$tpr)
    # part_auroc_normed <- trapz(min_max_normalization(roc_tbl_fpr_1e3$fpr), 
    #                     min_max_normalization(roc_tbl_fpr_1e3$tpr))
    # use spline smoothing to find the TPR value at FPR == 1e-3
    # ROC_splinefun <- splinefun(x = ROC_res_df$fpr, ROC_res_df$tpr, ties = "max")
    # tpr_at_fpr_1e_3 <- ROC_splinefun(1e-3)
    for (t in c("1e-2", "1e-3", "1e-4", "1e-5", "1e-6")){
      assign(x = paste0("pauroc_", gsub("-", "", t)),
             value = 
               spline_merge_pAUROC(ROC_res_df, ROC_spline_df, 
                                   as.numeric(t)))
    }
    # calculate pAUROC at different FPR values
    pauroc_df <- tibble(ROC_fil = roc_path,
                 pauroc_1e2 = pauroc_1e2,
                 pauroc_1e3 = pauroc_1e3,
                 # partial_auroc_1e3_norm = part_auroc_normed,
                 # tpr_at_fpr_1e_3 = tpr_at_fpr_1e_3,
                 pauroc_1e4 = pauroc_1e4,
                 pauroc_1e5 = pauroc_1e5,
                 pauroc_1e6 = pauroc_1e6,
                 # pauroc_1e7 = pauroc_1e7,
                 # pauroc_1e8 = pauroc_1e8,
                 # pauroc_1e9 = pauroc_1e9
    )
    write_csv(pauroc_df, pauroc_df_path)
  } else{
    pauroc_df <- read_csv(pauroc_df_path)
    pauroc_df
  }
}

## import results summary, including Arttu's calls
selex_results <- read_csv(file.path(proj_dir, "/../TEHMM_HT-SELEX_Results.csv")) %>%
  # mutate(selex_plate = gsub("(YW[TUV])_.*", "\\1", Exp_cycle1_identifier)) %>%
  # filter(Exp_Type == "HT-SELEX") %>%
  filter(!if_all(everything(), is.na)) %>%
  mutate(Arttu_call = ifelse(is.na(Arttu_call) & Exp_Type == "HT-SELEX", "fail", 
                             Arttu_call)) # assume NAs (uneval, altho I think he did look at them) are fails

## join with selex_results
selex_res_bmark <- left_join(selex_results, benchmark_summ_df, 
                             by = c("pTH", "Well", "selex_plate")) %>% 
  mutate(ROC_fil = paste0(exp_id, "/", motif, "_ROC.tsv"),
         Exp_Id = ifelse(is.na(Exp_Id), paste0(pTH, "_", selex_plate), Exp_Id))

## calculate partial AUROC values
ROC_fils <- selex_res_bmark %>% filter(!is.na(exp_id)) %>% pull(ROC_fil) %>%
  unique()

plan(multicore)
part_auroc_df <- ROC_fils[range_start:range_end] %>% future_map(part_auroc) %>% list_rbind()
plan(sequential)

selex_res_bmark_part_auroc <- right_join(selex_res_bmark, part_auroc_df) %>%
  mutate(Arttu_call = ifelse(is.na(Arttu_call), "fail", Arttu_call))

## save to file
write_csv(selex_res_bmark_part_auroc, 
          file.path(proj_dir, paste0("../selex_res_motif_bmark_", 
                    range_start, "_", range_end, ".csv.gz")))
