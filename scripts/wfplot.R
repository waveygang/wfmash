#!/usr/bin/env Rscript

library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
path_input_tsv <- args[1]

xxx <- read.table(path_input_tsv, sep = '\t', header = T)
xxx$info = as_factor(xxx$info)

p <- ggplot(xxx, aes(x=h, y=v,fill=info)) +
  geom_tile() +
  coord_fixed(ratio=1) + ggtitle("High-order DP matrix") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  #+ xlim(0,5000) + ylim(0, 5000)


path_output_png <- paste0(path_input_tsv, ".png")
ggsave(path_output_png, width = 11, height = 11, dpi = 600, units = "in", device='png', limitsize = FALSE)
