library(tidyverse)

# Convert Bed 3P-Seq files into point features for Stringtie2

args = commandArgs(trailingOnly=TRUE)

cpas <- read_tsv(args[1], col_names=c("chr","start","end","Name", "Score", "strand")) %>%
		mutate(pos = end)

cpas$pos[cpas$strand == "-"] <- cpas$start[cpas$strand == "-"]

cpas %>% mutate(cp = "CPAS") %>%
	dplyr::select(chr, pos, strand, cp) %>%
	write_tsv(file=gsub(".bed",".ptf", arg[1]), col_names=F)
