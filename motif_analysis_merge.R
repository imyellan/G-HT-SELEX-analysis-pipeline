library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(furrr)
library(pracma)
library(magrittr)
library(tidyr)

proj_dir <<- ifelse(Sys.info()["sysname"] == "Darwin",
                    "~/iyellan_bc2-clust/TEHMM_proj/analysis/selex_results/motif_pipeline",
                    "/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline")

bmark_pauroc_dfs <- list.files(
  file.path(proj_dir, ".."), pattern = "selex_res_motif_bmark_", full.names = T)
selex_res_bmark_part_auroc <- 
  bmark_pauroc_dfs %>% map(~read_csv(., col_types = c("Arttu_comment" = "c"))) %>% 
  list_rbind() %>% unique() %>% pivot_wider(names_from = "metric", values_from = "score") 

## annotate with TEHMM data
controls <- c("turtle_ZBED6" = "ENSPSIG00000001968", 
              "platypus_TIGD3" = "ENSOANG00000049052",
              "platypus_FLYWCH" = "ENSOANG00000028773",
              "croc_THAP1" = "ENSCPRG00005014599",
              "brugia_malayi_nematode_pax5_like" = "WBGene00222217",
              "zebrafish_naif1" = "ENSDARG00000062753",
              "elephant_shark_rag1" = "ENSCMIG00000015975"
) 
TE_selection <- read_csv("~/iyellan_bc2-clust/TEHMM_proj/analysis/proteome_scan/ensembl/TE_selection_seqs_order.csv") %>% transmute(
  genid = Name,
  exp_role = "TE",
  common = "TE",
  protid = gsub("_.*", "", Name),
  TE_fam = Name,
  `Descriptive name` = `Name`
)
all_TEHMM_results_df <- read_csv("~/iyellan_bc2-clust/TEHMM_proj/analysis/proteome_scan/ensembl/all_scan_df_pfam_passed_ogs_dom_TE_sets_sorted_treedists_sel_af_intrus.csv.gz")
final_selection_df <- read_csv("~/iyellan_bc2-clust/TEHMM_proj/analysis/proteome_scan/ensembl/final_selection_flanks_v5.csv") %>% select(all_of(1:17)) %>% 
  left_join(all_TEHMM_results_df %>% select(protid, genid, TE_fam)) %>%
  mutate(
    exp_role = ifelse(genid %in% controls, "control", NA),
    ensembl_gene_name_up = ifelse(is.na(ensembl_gene_name_up), orthogroup_nm, ensembl_gene_name_up),
    `Descriptive name` = 
           gsub(" ", "_", paste(protid, ensembl_gene_name_up, TE_fam, common, sep = "_"))) %>%
  select(-c(unip_length, comb_length, env_from, env_to, ens_length.x, ens_length.y, uniprot,
            uniprot_nm)) %>% unique() 
final_selection_df_TEs <- bind_rows(final_selection_df, TE_selection)
selex_res_bmark_part_auroc <- selex_res_bmark_part_auroc %>% 
  left_join(final_selection_df_TEs) %>%
  mutate(exp_role = case_when(
    grepl("ZBED6|ZBED8|ZNF862|THAP3|HARBI1|NAIF1|MSATND2|TSNARE1|PGBD4|SETMAR|TIGD6|JRKL|TIGD2|POGZ|ZMYM1|ZMYM5|ZMYM6",
          ensembl_gene_name_up) ~ "previously_IDd_human",
    is.na(exp_role) ~ "novel",
    .default = exp_role))

## get top motif score, plot comparison to Arttu's calls
# score_summ_p <-
#   selex_res_bmark_part_auroc %>% 
#   group_by(Exp_Id) %>% 
#   tidyr::pivot_wider(names_from = "metric", values_from = "score") %>%
#   filter(!is.na(ROC), !is.na(PR)) %>%
#   slice_max(tibble(ROC, PR), n = 1, with_ties = F) %>%
#   ggplot(aes(x = ROC, y = partial_auroc_norm)) + 
#   geom_bin2d() + 
#   facet_wrap(vars(Arttu_call))
# ggsave(score_summ_p, file = file.path(proj_dir, "selex_res_motif_bmark_top.pdf"))

## same but using the partial roc score
for (t in 2:6) {
  pauroc_thresh <-
    case_when(t == 4 ~ 1e-6, 
              t == 5 ~ 1e-7, 
              t == 2 ~ 0.000316,
              t == 6 ~1e-8,
              .default = 1e-5)
  top_motifs_df <-
    selex_res_bmark_part_auroc %>%
    mutate(#tpr_at_fpr_1e_3 = ifelse(tpr_at_fpr_1e_3 > 1, 1, tpr_at_fpr_1e_3),
      Arttu_call = paste0("Arttu: ", Arttu_call),
      # high_part_auroc = ifelse(partial_auroc_norm >= 0.85, "pAUROC >= 0.85", "pAUROC < 0.85")
    ) %>%
    group_by(Exp_Id) %>% 
    filter(!is.na(ROC), !is.na(PR)) %>%
    mutate(
      good_ROC = ROC >= 0.85,
      top_motif = 
        (good_ROC & ROC == max(ROC)) | 
        (!good_ROC & !!sym(paste0("pauroc_1e",t)) == max(!!sym(paste0("pauroc_1e",t)))),
    ) %>%
    filter(top_motif) %>% 
    slice_max(PR) %>% slice_sample(n = 1) %>% # for handling ties
    mutate(log10_pauroc = log10(!!sym(paste0("pauroc_1e",t))))
  label_df <- top_motifs_df %>% 
    filter((Arttu_call == "Arttu: ok" & !!sym(paste0("pauroc_1e",t)) < pauroc_thresh) |
             ((Arttu_call == "Arttu: fail" | Arttu_call == "Arttu: maybe") & 
                !!sym(paste0("pauroc_1e",t)) > pauroc_thresh))
    top_motifs_df %>%
    mutate(log10_pauroc = log10(!!sym(paste0("pauroc_1e",t)))) %>%
      ggplot(aes(x = log10_pauroc, y = ROC, colour = Arttu_call)) + 
      geom_count(alpha = 0.7) +  
      scale_size_area(max_size = 7) +
      geom_label_repel(data = label_df %>% mutate(short_lab = substr(`Descriptive name`, 1, 15)), 
                      aes(x = log10_pauroc, y = ROC, label = short_lab),
                      min.segment.length = 0.8) +
    # facet_wrap(vars(high_part_auroc)) + 
    theme_bw() +
    geom_hline(yintercept = 0.85, linetype = "dashed") +
    geom_vline(xintercept = log10(pauroc_thresh), linetype = "dashed", colour = "red") +
    labs(x = paste0("Log10 pAUROC at FPR = 1E-", t), y = "AUROC")# +
  # coord_cartesian(xlim = c(-10, -4))
  ggsave(file = file.path(proj_dir, paste0("selex_res_motif_bmark_top_pauroc_1e-", t, ".pdf")),
         width = 8, height = 6)
}

## save df of working and non-working motifs
top_motifs_dir <- file.path(proj_dir, "../top_working_motifs")
dir.create(top_motifs_dir, showWarnings = F, recursive = T)

top_motifs_df <-
  selex_res_bmark_part_auroc %>%
  mutate(Arttu_call = paste0("Arttu: ", Arttu_call)) %>%
  group_by(Exp_Id) %>% 
  filter(!is.na(ROC), !is.na(PR)) %>%
  mutate(
    good_ROC = ROC >= 0.85,
    top_motif = 
      (good_ROC & ROC == max(ROC)) | 
      (!good_ROC & pauroc_1e4 == max(pauroc_1e4)),
  ) %>%
  filter(top_motif) %>%
  slice_max(PR)

top_working_motifs_df <- top_motifs_df %>% 
  mutate(working_experiment = good_ROC | (!good_ROC & pauroc_1e4 >= 1E-6) | 
           `Descriptive name` == "hAT-31_DR_tp#DNA/hAT-Ac" |
           `Descriptive name` == "ENSIPUP00000009670.1_PAX1B_TcMar-ISRm11_Channel_catfish"
           )
sort_working_dfs <- function(df){
  out_dir <- 
    gsub(" |\\(|\\)|#|\\||;|'|:", "_", file.path(top_motifs_dir, 
                                        gsub("/", "_", paste0(df$`Descriptive name`, "_", 
                                                    df$Exp_Type, "_", df$selex_plate))))
  dir.create(out_dir, showWarnings = F)
  
  base_dir <- file.path(proj_dir, dirname(df$ROC_fil))
  motif_ROC_figs <- list.files(base_dir, pattern = basename(gsub(".tsv", ".*.pdf", df$ROC_fil)), 
                               full.names = T)
  motif_logo_fils <- list.files(base_dir, 
                                pattern = gsub("ROC.*", "logo.*.png", basename(motif_ROC_figs)),
                                full.names = T)
  autoseed_base_dir <- ifelse(df$selex_plate == "YWT", 
                              "/Users/isaacyellan/rc_home-clust/hughespub/finishTF/YWT_Mostly_Isaac/AutoseedOutput",
                              "/Users/isaacyellan/rc_home-clust/hughespub/SELEX_Data/Autoseed_outputs")
  autoseed_fils <- list.files(autoseed_base_dir, pattern = paste0(".*", df$pTH, ".*.svg"),
                              full.names = T)
  file.copy(c(motif_ROC_figs, motif_logo_fils, autoseed_fils), paste0(out_dir, "/"), copy.date = T)
}
top_working_motifs_df %>% 
  rowwise() %>%
  group_split() %>%
  walk(sort_working_dfs)

