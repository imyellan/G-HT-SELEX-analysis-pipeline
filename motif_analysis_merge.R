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
  "/home/hugheslab1/iyellan/TEHMM_proj/analysis/selex_results/motif_pipeline"
)


sort_working_dfs <- function(df) {
  top_motif <- ifelse(df$top_motif, "TOP_", "")
  out_dir <-
    gsub(
      " |\\(|\\)|#|\\||;|'|:", "_",
      file.path(
        top_motifs_dir,
        gsub("/", "_", paste0(
          df$`Descriptive name`, "_",
          df$Exp_Type, "_", df$selex_plate
        ))
      )
    )
  dir.create(out_dir, showWarnings = F)

  base_dir <- file.path(proj_dir, dirname(df$ROC_fil_0.1))
  motif_ROC_figs <-
    list.files(base_dir,
      pattern = gsub(
        "[0-9]+_ROC", ".*_ROC",
        basename(gsub(".tsv", ".pdf", df$ROC_fil_0.1))
      ),
      full.names = T
    )
  motif_ROC_figs <- grep("partial", motif_ROC_figs, value = T, invert = T)
  motif_logo_fils <- list.files(base_dir,
    pattern = unique(gsub(
      "[0-9]+_ROC.*", ".*logo.*.png",
      basename(df$ROC_fil_0.1)
    )),
    full.names = T
  )
  motif_logo_fils <- grep(".png", motif_logo_fils, value = T)
  autoseed_pref <- ifelse(Sys.info()["sysname"] == "Darwin", 
                          "/Users/isaacyellan/rc_home-clust/hughespub",
                          "/home/hughespub")
  autoseed_base_dir <- ifelse(df$selex_plate == "YWT",
    file.path(autoseed_pref, "finishTF/YWT_Mostly_Isaac/AutoseedOutput"),
    file.path(autoseed_pref, "SELEX_Data/Autoseed_outputs")
  )
  autoseed_fils <- list.files(autoseed_base_dir,
    pattern = paste0(".*", df$pTH, ".*.svg"),
    full.names = T
  )
  file.copy(autoseed_fils, paste0(out_dir, "/"), copy.date = T)
  motif_pref <- dirname(df$ROC_fil_0.1)
  if(length(motif_logo_fils) >= 1){
    for (fil in motif_logo_fils) {
      file.copy(fil, paste0(out_dir, "/", top_motif, motif_pref, "_", basename(fil)))
    }
  }
  if(length(motif_ROC_figs) >= 1){
    for (fil in motif_ROC_figs) {
      file.copy(fil, paste0(out_dir, "/", top_motif, motif_pref, "_", basename(fil)))
    }
  }
}

bmark_pauroc_dfs <- list.files(
  file.path(proj_dir, ".."),
  pattern = "selex_res_motif_bmark_", full.names = T
)
selex_res_bmark_part_auroc <-
  bmark_pauroc_dfs %>%
  map(~ read_csv(., col_types = c("Arttu_comment" = "c"))) %>%
  list_rbind() %>%
  unique() # %>%
# pivot_wider(names_from = "metric", values_from = "score", )

## annotate with TEHMM data
controls <- c(
  "turtle_ZBED6" = "ENSPSIG00000001968",
  "platypus_TIGD3" = "ENSOANG00000049052",
  "platypus_FLYWCH" = "ENSOANG00000028773",
  "croc_THAP1" = "ENSCPRG00005014599",
  "brugia_malayi_nematode_pax5_like" = "WBGene00222217",
  "zebrafish_naif1" = "ENSDARG00000062753",
  "elephant_shark_rag1" = "ENSCMIG00000015975"
)
TE_selection <- read_csv(file.path(proj_dir, "../../proteome_scan/ensembl/TE_selection_seqs_order.csv")) %>% transmute(
  genid = Name,
  exp_role = "TE",
  common = "TE",
  protid = gsub("_.*", "", Name),
  TE_fam = Name,
  `Descriptive name` = `Name`
)
all_TEHMM_results_df <- read_csv(file.path(proj_dir, "../../proteome_scan/ensembl/all_scan_df_pfam_passed_ogs_dom_TE_sets_sorted_treedists_sel_af_intrus.csv.gz"))
final_selection_df <- read_csv(file.path(proj_dir, "../../proteome_scan/ensembl/final_selection_flanks_v5.csv")) %>%
  select(all_of(1:17)) %>%
  left_join(all_TEHMM_results_df %>% select(protid, genid, TE_fam, clade)) %>%
  mutate(
    exp_role = ifelse(genid %in% controls, "control", NA),
    ensembl_gene_name_up = ifelse(is.na(ensembl_gene_name_up), orthogroup_nm, ensembl_gene_name_up),
    `Descriptive name` =
      gsub(" ", "_", paste(protid, ensembl_gene_name_up, TE_fam, common, sep = "_"))
  ) %>%
  select(-c(
    unip_length, comb_length, env_from, env_to, ens_length.x, ens_length.y, uniprot,
    uniprot_nm
  )) %>%
  unique()
final_selection_df_TEs <- bind_rows(final_selection_df, TE_selection)
selex_res_bmark_part_auroc <- selex_res_bmark_part_auroc %>%
  left_join(final_selection_df_TEs) %>%
  mutate(
    ensembl_gene_name_up = case_when(
      pTH == "UT380-226" ~ "TIGD4",
      pTH == "UT380-242" ~ "ZBED2",
      pTH == "UT380-245" ~ "ZBED5",
      .default = ensembl_gene_name_up
    ),
    TE_fam = case_when(
      pTH == "UT380-226" ~ "TcMar-Tigger",
      pTH == "UT380-242" ~ "hAT-Ac",
      pTH == "UT380-245" ~ "hAT-Charlie",
      .default = TE_fam
    ),
    common = ifelse(grepl("UT380", pTH), "Human", common),
    exp_role = case_when(
      grepl(
        "ZBED6|ZBED8|ZNF862|THAP3|HARBI1|NAIF1|MSANTD2|TSNARE1|PGBD4|SETMAR|TIGD6|JRKL|TIGD2|POGZ|ZMYM1|ZMYM5|ZMYM6|GTF2IRD2|TIGD4|ZBED2|ZBED5",
        ensembl_gene_name_up
      ) & common == "Human" ~ "Known human TE-derived",
      is.na(exp_role) ~ "novel",
      .default = exp_role
    ),
    clade = ifelse(grepl("Cat|Rabbit|Human", common), "Mammals", 
                   ifelse(exp_role == "TE", "TE", clade)),
    input_type = ifelse(grepl("PEAK", exp_id), "MAGIX",
      ifelse(grepl("GHT", exp_id), "GHT reads",
        "HT reads"
      )
    ),
    genid = 
      case_when(
        protid_stab == "ENSFCAP00000042068" ~ "ENSFCAG00000040371",
        protid_stab == "ENSOCUP00000042417" ~ "ENSOCUG00000007931",
        `Descriptive name` == "TIGD4-DBD" ~ "ENSG00000169989",
        `Descriptive name` == "ZBED2-FL" ~ "ENSG00000177494",
        `Descriptive name` == "ZBED5-DBD" ~ "ENSG00000236287",
        .default = genid
      )
  )

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

selex_res_bmark_part_auroc %>%
  filter(grepl("C17ORF113", `Descriptive name`)) %>%
  filter(metric == "PR") %>%
  mutate(
    exp_motif = paste(exp_id, motif, sep = "_"),
    motif_caller = gsub("_.*", "", motif)
  ) %>%
  ggplot(aes(
    x = bmark_pos_num, y = score, group = exp_motif,
    colour = motif_caller
  )) +
  geom_line(show.legend = T) +
  facet_wrap(vars(input_type)) +
  scale_x_log10()

## NEW METHOD FOR SELECTING THAT JUST USES THE ROCS AT DIFFERENT POSITIVE FRACS
motif_eval_df <-
  selex_res_bmark_part_auroc %>%
  filter(input_type != "GHT reads") %>% # something funky about these
  group_by(exp_id, motif) %>%
  # motif must be benchmarked at all 3 positive fractions, & each fraction must be >= 10 sequences (reads or peaks) for benchmarking
  mutate(
    # enough_seqs =
    #   n_distinct(bmark_pos_frac) == 3 & 
    #   (all(bmark_pos_num >= 10) | any(input_type == "MAGIX" & bmark_pos_num >= 5 & bmark_pos_frac == 0.1)),
    Arttu_call = paste0("Arttu: ", Arttu_call)
  ) %>%
  group_by(Exp_Id) %>%
  filter(all(!is.na(metric))) %>%
  select(-contains("pauroc")) %>%
  pivot_wider(names_from = metric, values_from = score) %>%
  pivot_wider(
    names_from = bmark_pos_frac,
    values_from = c(ROC, PR, bmark_pos_num, ROC_fil)
  ) %>%
  mutate(
    across(contains("ROC_0"), ~ . > 0.8, .names = "{.col}_good_ROC"),
    across(contains("PR_0"), ~ . > 0.8, .names = "{.col}_good_PR")
  ) %>%
  mutate(good_motif = 
           if_all(matches("ROC_0\\.[015]"), ~!is.na(.)) &
           ((ROC_0.1_good_ROC & bmark_pos_num_0.1 >= 10) | 
           ((ROC_0.01_good_ROC & PR_0.01_good_PR) & bmark_pos_num_0.01 >= 5))) %>%
  group_by(`Descriptive name`, input_type) %>%
  mutate(
    good_ROC = ifelse(good_motif, ROC_0.1, NA)) %>%
  arrange(`Descriptive name`, input_type, desc(good_ROC), desc(PR_0.1), 
          desc(ROC_0.5), desc(PR_0.5)) %>%
  mutate(top_motif = good_motif & (row_number() == 1))
write_csv(motif_eval_df, file.path(proj_dir, "../motif_evaluation_summ_metadat.csv.gz"))

date_str <- system("date +%y%m%d", intern = T)
top_motifs_dir <- file.path(proj_dir, "../top_working_motifs", date_str)
dir.create(top_motifs_dir, showWarnings = F, recursive = T)
# system(paste0("rm -r ", top_motifs_dir, "/*"))
motif_eval_df %>%
  filter(good_motif) %>%
  # group_by(`Descriptive name`, Exp_Type) %>%
  rowwise() %>%
  group_split() %>%
  walk(sort_working_dfs)

## same but using the partial roc score
# for (t in 2:6) {
#   pauroc_thresh <-
#     case_when(t == 4 ~ 1e-6,
#               t == 5 ~ 1e-7,
#               t == 2 ~ 0.000316,
#               t == 6 ~1e-8,
#               .default = 1e-5)
#
#       good_ROC = ROC >= 0.85,
#       top_motif =
#         (good_ROC & ROC == max(ROC)) |
#         (!good_ROC & !!sym(paste0("pauroc_1e",t)) == max(!!sym(paste0("pauroc_1e",t)))),
#     ) %>%
#     filter(top_motif) %>%
#     slice_max(PR) %>% slice_sample(n = 1) %>% # for handling ties
#     mutate(log10_pauroc = log10(!!sym(paste0("pauroc_1e",t))))
#   label_df <- top_motifs_df %>%
#     filter((Arttu_call == "Arttu: ok" & !!sym(paste0("pauroc_1e",t)) < pauroc_thresh) |
#              ((Arttu_call == "Arttu: fail" | Arttu_call == "Arttu: maybe") &
#                 !!sym(paste0("pauroc_1e",t)) > pauroc_thresh))
#     top_motifs_df %>%
#     mutate(log10_pauroc = log10(!!sym(paste0("pauroc_1e",t)))) %>%
#       ggplot(aes(x = log10_pauroc, y = ROC, colour = Arttu_call)) +
#       geom_count(alpha = 0.7) +
#       scale_size_area(max_size = 7) +
#       geom_label_repel(data = label_df %>% mutate(short_lab = substr(`Descriptive name`, 1, 15)),
#                       aes(x = log10_pauroc, y = ROC, label = short_lab),
#                       min.segment.length = 0.8) +
#     # facet_wrap(vars(high_part_auroc)) +
#     theme_bw() +
#     geom_hline(yintercept = 0.85, linetype = "dashed") +
#     geom_vline(xintercept = log10(pauroc_thresh), linetype = "dashed", colour = "red") +
#     labs(x = paste0("Log10 pAUROC at FPR = 1E-", t), y = "AUROC")# +
#   # coord_cartesian(xlim = c(-10, -4))
#   ggsave(file = file.path(proj_dir, paste0("selex_res_motif_bmark_top_pauroc_1e-", t, ".pdf")),
#          width = 8, height = 6)
# }
#
# ## save df of working and non-working motifs
# date_str <- system("date +%y%m%d", intern = T)
# top_motifs_dir <- file.path(proj_dir, "../top_working_motifs", date_str)
# dir.create(top_motifs_dir, showWarnings = F, recursive = T)
#
# top_motifs_df <-
#   selex_res_bmark_part_auroc %>%
#   mutate(Arttu_call = paste0("Arttu: ", Arttu_call)) %>%
#   group_by(Exp_Id) %>%
#   filter(!is.na(ROC), !is.na(PR)) %>%
#   mutate(
#     good_ROC = ROC >= 0.85,
#     top_motif =
#       (good_ROC & ROC == max(ROC)) |
#       (!good_ROC & pauroc_1e4 == max(pauroc_1e4)),
#   ) %>%
#   filter(top_motif) %>%
#   slice_max(PR)
# top_working_motifs_df <- top_motifs_df %>%
#   mutate(working_experiment = good_ROC | (!good_ROC & pauroc_1e4 >= 1E-6) |
#            `Descriptive name` == "hAT-31_DR_tp#DNA/hAT-Ac" |
#            `Descriptive name` == "ENSIPUP00000009670.1_PAX1B_TcMar-ISRm11_Channel_catfish"
#   )

#
