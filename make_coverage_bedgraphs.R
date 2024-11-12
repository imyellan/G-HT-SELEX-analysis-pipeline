library(dplyr)
library(tidyr)
library(readr)
library(dtplyr)

coverage_df <- read_tsv(commandArgs(trailingOnly = T)[1])
outdir <- commandArgs(trailingOnly = T)[2]
gene <- basename(outdir)

coverage_df <- coverage_df %>%
  lazy_dt() %>%
  pivot_longer(-all_of(1:14), names_to = "cycle_num_pos", values_to = "coverage") %>%
  mutate(
    type = gsub("_cycle.*", "", cycle_num_pos),
    cycle_num = gsub(".*_(cycle[0-9])_.*", "\\1", cycle_num_pos),
    rel_pos = as.numeric(gsub(".*_([0-9]+)$", "\\1", cycle_num_pos)),
    new_start = start + rel_pos - 1,
    new_end = start + rel_pos
  ) %>%
  select(X.chr, new_start, new_end, coverage, type, cycle_num) %>%
  as_tibble()

coverage_df %>%
  group_by(type, cycle_num) %>%
  group_walk(~ .x %>%
    arrange(X.chr, new_start) %>%
    write_tsv(
      file =
        file.path(
          outdir,
          paste0(
            gene, "_coverage_",
            .y$type, "_", .y$cycle_num, ".bed"
          )
        ),
      col_names = F
    ))
