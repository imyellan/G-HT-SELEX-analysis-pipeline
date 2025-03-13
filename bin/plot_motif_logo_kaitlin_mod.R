#!/usr/bin/env Rscript

library(ggseqlogo)
library(ggplot2)
library(universalmotif)

args = commandArgs(TRUE)
filename = args[1]
outdir = "."
# outprefix = paste0(gsub('.txt','', basename(filename)), "_logo")
outprefix <- "logo"

# pfm <- read.table(filename, sep = "\t", stringsAsFactors = FALSE,
#                   header = FALSE, skip = 3,row.names = 1)
pfm <- read_matrix(filename, positions = "rows", skip = 0, headers = F,
                   alphabet = "DNA")@motif
#pfm <- t(pfm)
#colnames(pfm) <- c("A", "C", "G", "T")
rownames(pfm) <- c("A", "C", "G", "T")
#pfm <- t(pfm)

colscheme <- make_col_scheme(chars = c("A", "C", "G", "T"), 
                             cols=c("#00CC00", "#0000CC", "#FFB302", "#CC0001"),
                             name="DNAalph")
plot <- ggplot() + 
  geom_logo(as.matrix(pfm), 
            namespace = colscheme$letter, 
            col_scheme = colscheme, 
            font="roboto_bold",
            method="bits") + 
  ylim(c(0, 2)) + 
  theme_logo() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(0, "null"),
        axis.ticks.margin = unit(0, "null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,2)) 

png(file.path(outdir, paste0(outprefix, ".png")), width=(1.5*ncol(pfm)/10), height=0.75, units='in', res=300)
print(plot)
dev.off()

rev_pfm <- pfm[, ncol(pfm):1]
colnames(rev_pfm) <- 1:ncol(rev_pfm)
rev_pfm <- rev_pfm[c("T", "G", "C", "A"), ]
rownames(rev_pfm) <- c("A", "C", "G", "T")

plot <- ggplot() + 
  geom_logo(as.matrix(rev_pfm), 
            namespace = colscheme$letter, 
            col_scheme = colscheme, 
            font="roboto_bold",
            method="bits") + 
  ylim(c(0, 2)) + 
  theme_logo() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(0, "null"),
        axis.ticks.margin = unit(0, "null")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,2)) 

png(file.path(outdir, paste0(outprefix, "_REV.png")), 
    width=(1.5*ncol(pfm)/10), height=0.75, units='in', res=300)
print(plot)
dev.off()
