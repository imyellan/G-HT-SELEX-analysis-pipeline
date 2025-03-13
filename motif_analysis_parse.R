#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(furrr)
library(pracma)
library(magrittr)

proj_dir <<- ifelse(Sys.info()["sysname"] == "Darwin",
  "~/iyellan_bc2-clust/TEHMM_proj/analysis/selex_results/motif_pipeline",
  "/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline"
)
range_start <- as.numeric(commandArgs(trailingOnly = T)[1])
range_end <- as.numeric(commandArgs(trailingOnly = T)[2])
recalc <- as.logical(commandArgs(trailingOnly = T)[3])
recalc <- ifelse(is.na(recalc), FALSE, recalc)
benchmark_summ_df <- read_csv(file.path(proj_dir, "../motif_benchmark_summary_df.csv.gz"))

import_ROC <- function(f) {
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

spline_merge_pAUROC <- function(ROC_res_df, ROC_spline_df, thresh) {
  ROC_res_df_thresh <- ROC_res_df %>% filter(fpr <= thresh)
  if (max(ROC_res_df_thresh$fpr) < thresh) {
    if (nrow(ROC_res_df_thresh) == 0) {
      ROC_res_df_thresh <- data.frame(fpr = 0, tpr = 0)
    }
    ROC_spline_df %<>% filter(fpr <= thresh)
    if (nrow(ROC_spline_df) == 0) {
      ROC_spline_df <- data.frame(fpr = 0, tpr = 0)
    }
    unique(rbind(
      ROC_res_df_thresh %>% select(fpr, tpr),
      ROC_spline_df %>% select(fpr, tpr)
    )) %$%
      trapz(fpr, tpr)
  } else {
    trapz(ROC_res_df_thresh$fpr, ROC_res_df_thresh$tpr)
  }
}

part_auroc <- function(roc_path, recalc = F) {
  pauroc_df_path <- gsub("_ROC.tsv", "_pauroc_df.csv.gz", file.path(proj_dir, roc_path))
  if (recalc #| !file.exists(pauroc_df_path)
  ) {
    ROC_res_df <- import_ROC(file.path(proj_dir, roc_path))
    ## while at it, save an ROC plot
    ROC_res_df %>%
      ggplot(aes(fpr, tpr)) +
      geom_line() +
      ggtitle(roc_path) +
      theme_bw() +
      theme(legend.position = "none")
    ggsave(gsub(".tsv", ".pdf", file.path(proj_dir, roc_path)))

    ## create a spline of the ROC curve
    ROC_spline <- spline(
      x = ROC_res_df$fpr, ROC_res_df$tpr,
      n = 1000 * length(ROC_res_df$fpr), ties = "max"
    )
    ROC_spline_df <- data.frame(fpr = ROC_spline$x, tpr = ROC_spline$y)

    ## save partial ROC plot - using FPR <= 1e-3, TPR <= 1e-2
    ROC_res_df %>%
      ggplot(aes(fpr, tpr)) +
      geom_line() +
      ggtitle(roc_path) +
      theme_minimal() +
      theme(legend.position = "none") +
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
    for (t in c("1e-2", "1e-3", "1e-4", "1e-5", "1e-6")) {
      assign(
        x = paste0("pauroc_", gsub("-", "", t)),
        value =
          spline_merge_pAUROC(
            ROC_res_df, ROC_spline_df,
            as.numeric(t)
          )
      )
    }
    # calculate pAUROC at different FPR values
    pauroc_df <- tibble(
      ROC_fil = roc_path,
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
  } else {
    if (file.exists(pauroc_df_path)) {
      pauroc_df <- read_csv(pauroc_df_path)
      pauroc_df
    } else {
      data.frame(
        ROC_fil = roc_path,
        pauroc_1e2 = NA,
        pauroc_1e3 = NA,
        pauroc_1e4 = NA,
        pauroc_1e5 = NA,
        pauroc_1e6 = NA
      )
    }
  }
}

## import results summary, including Arttu's calls
selex_results <- read_csv(file.path(proj_dir, "/../TEHMM_HT-SELEX_Results.csv")) %>%
  # mutate(selex_plate = gsub("(YW[TUV])_.*", "\\1", Exp_cycle1_identifier)) %>%
  # filter(Exp_Type == "HT-SELEX") %>%
  filter(!if_all(everything(), is.na)) %>%
  mutate(Arttu_call = ifelse(is.na(Arttu_call) & Exp_Type == "HT-SELEX", "fail",
    Arttu_call
  )) # assume NAs (uneval, altho I think he did look at them) are fails

## deal with MAGIX peaks-based motifs
match_genes <- function(gene_nm) {
  match <- grep(gene_nm, unique(selex_results$`Descriptive name`), value = T)
  if (length(match) == 1) {
    tibble(gene = gene_nm, `Descriptive name` = match)
  } else {
    tibble(gene = gene_nm, `Descriptive name` = NA) %>%
      mutate(
        `Descriptive name` =
          case_when(
            gene == "ACOM031257_PA_1" ~
              "ACOM031257-PA.1_151571_4_hAT-Tip100_Anopheles_coluzzii",
            gene == "NAIF1" ~
              "ENSP00000362170.4_NAIF1_PIF-Harbinger_Human",
            gene == "POGZ" ~
              "ENSP00000357856.2_POGZ_TcMar-Tc2_Human",
            gene == "UTF1" ~
              "ENSP00000305906.2_UTF1_PIF-Harbinger_Human",
            gene == "ZBED6" ~
              "ENSP00000447879.1_ZBED6_hAT-Ac_Human",
            gene == "ZNF862" ~
              "ENSP00000223210.4_ZNF862_hAT-Tip100_Human"
          )
      )
  }
}
peak_rocs <- benchmark_summ_df %>%
  filter(cycle == "MAGIX") %>%
  mutate(gene = gsub("target_|_LTR.*", "", pTH))
peak_rocs <- peak_rocs %>%
  left_join(peak_rocs %>% pull(gene) %>% unique() %>%
    map(match_genes) %>% list_rbind()) %>%
  select(-gene)
selex_res_bmark_peaks <-
  left_join(
    selex_results %>%
      filter(Exp_Type == "GHT-SELEX", selex_plate != "YWT") %>%
      select(-Well, -Exp_cycle1_identifier, -Exp_Id) %>% unique(),
    peak_rocs %>% select(-pTH, -selex_plate) %>% unique()
  )

## join with selex_results
selex_res_bmark_reads <- left_join(selex_results, benchmark_summ_df %>%
  filter(cycle != "MAGIX"),
by = c("pTH", "Well", "selex_plate")
)
selex_res_bmark <- bind_rows(selex_res_bmark_reads, selex_res_bmark_peaks) %>%
  unique() %>%
  mutate(
    ROC_fil = paste0(
      exp_id, "/", motif, "_",
      gsub("\\.", "", as.character(bmark_pos_frac)), "_ROC.tsv"
    ),
    Exp_Id = ifelse(is.na(Exp_Id), paste0(pTH, "_", selex_plate), Exp_Id)
  )

## calculate partial AUROC values
ROC_fils <- selex_res_bmark %>%
  filter(!is.na(exp_id)) %>%
  pull(ROC_fil) %>%
  unique()

plan(multicore)
part_auroc_df <- ROC_fils[range_start:range_end] %>%
  future_map(~ part_auroc(.x, recalc)) %>%
  list_rbind()
plan(sequential)

selex_res_bmark_part_auroc <- right_join(selex_res_bmark, part_auroc_df) %>%
  mutate(Arttu_call = ifelse(is.na(Arttu_call), "fail", Arttu_call))

## save to file
write_csv(
  selex_res_bmark_part_auroc,
  file.path(proj_dir, paste0(
    "../selex_res_motif_bmark_",
    range_start, "_", range_end, ".csv.gz"
  ))
)
