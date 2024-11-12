library(readr)
library(tidyr)
library(ggplot2)

proj_dir <- ifelse(Sys.info()["sysname"] == "Darwin",
                   "~/iyellan_bc2-clust/TEHMM_proj/analysis/selex_results/motif_pipeline",
                   "/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline")

bad_homer_roc <- read_tsv(file.path(proj_dir, "YWU_A_1_C08_pTH14640_CC40NCCAATA_eGFP_IVT_S50_R1_001/homer_out_10_2_ROC.tsv")) %>%
  mutate(pwm = "bad_homer")
good_streme_roc <- read_tsv(file.path(proj_dir, "YWU_A_3_C08_pTH14640_CC40NCCAATA_eGFP_IVT_S434_R1_001/streme_pfm_2_ROC.tsv")) %>%
  mutate(pwm = "good_streme")

roc_comb <- rbind(bad_homer_roc, good_streme_roc)
roc_comb %>% 
  ggplot(aes(x = fpr, y = tpr, color = pwm)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "PWM") +
  ggtitle("ROC curve for HOMER and STREME PWMs")

auc(bad_homer_roc$tpr, bad_homer_roc$fpr)


### failers
plot_ROCs <- function(in_tbl, key){
  pTH <- unique(in_tbl$pTH)
  import_ROC <- function(f){
    read_tsv(f) %>% 
      mutate(base = gsub("_ROC.tsv", "", gsub(".*motif_pipeline/", "", f)))
  }
  ROC_res <- list.files(
    list.files(proj_dir, pattern = pTH, include.dirs = T, full.names = T),
    pattern = "ROC.tsv", full.names = T, recursive = T)
  ROC_res_df <- ROC_res %>% map(import_ROC) %>% list_rbind()
  partial_int <- function(roc_tbl, key){
    min_max_normalization <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))
    }
    # browser()
    roc_tbl_low_fpr <- roc_tbl %>% filter(fpr <= 1e-3)
    tibble(partial_int_norm = trapz(min_max_normalization(roc_tbl_low_fpr$fpr), 
                               min_max_normalization(roc_tbl_low_fpr$tpr)))
  }
  
  ROC_res_part_int_df <- 
    ROC_res_df %>% group_by(base) %>% group_modify(~partial_int(.x,.y))
  high_part_int <- ROC_res_part_int_df %>% filter(partial_int_norm >= 0.85) %>% pull(base)
  # top tpr at low fpr
  p_high_tpr <- ROC_res_df %>% 
    # group_by(base) %>% 
    filter(base %in% high_part_int) %>%
    # filter(grepl("streme", base)) %>%
    ggplot(aes(x = fpr, y = tpr, colour = base)) +
    geom_line(linewidth = 0.1) +
    theme_bw() +
    guides(colour = guide_legend(nrow = 3)) +
    theme(legend.position = "top", legend.text = element_text(size = 7))
  p_low_tpr <- ROC_res_df %>% 
    group_by(base) %>% 
    filter(!base %in% high_part_int) %>%
    # filter(grepl("streme", base)) %>%
    ggplot(aes(x = fpr, y = tpr, colour = base)) +
    geom_line(linewidth = 0.1, show.legend = F) +
    theme_bw()
  p_high_tpr/p_low_tpr 
  ggsave(file.path("~/Downloads", paste0(pTH, ".pdf")), width = 10)
}

ok_ROC_mismatch <- selex_res_bmark %>% 
  group_by(Exp_Id) %>% 
  tidyr::pivot_wider(names_from = "metric", values_from = "score") %>%
  filter(all(ROC <= 0.8), Arttu_call == "ok")

ok_ROC_mismatch %>%
  group_walk(~plot_ROCs(.x,.y))


