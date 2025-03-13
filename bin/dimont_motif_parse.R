#!/usr/bin/env Rscript

library(readr)
library(xml2)

in_xml <- commandArgs(trailingOnly = T)[1]
pwm_num <- gsub(".*-pwm-([0-9]+).xml", "\\1", in_xml)

xml_obj <- read_xml(in_xml)
xml_list <- as_list(xml_obj)
nodeset <- xml_children(xml_obj)
motif_length <- xml_double(xml_contents(nodeset[2]))

pwm_df <- data.frame()
for (i in 3:(motif_length+2)) {
  pwm_df <- rbind(pwm_df, xml_double(xml_contents((xml_children(nodeset[i]))))[3:6])
}

write_tsv(pwm_df, file.path(dirname(in_xml), 
                            paste0("dimont_ght_pfm_", pwm_num, ".pfm")),
                            col_names = F)
