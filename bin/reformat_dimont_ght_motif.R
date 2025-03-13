#!/usr/bin/env Rscript

library(universalmotif)
library(readr)

motif_path <- commandArgs(trailingOnly = T)[1]

in_motif <- read_matrix(motif_path, type = "PWM", positions = "columns", rownames = T)
in_motif_ppm <- convert_type(in_motif, type = "PPM")

## save motif
write_tsv(as.data.frame(t(in_motif_ppm@motif)), gsub("pwm", "pfm", motif_path),
    col_names = F
)
