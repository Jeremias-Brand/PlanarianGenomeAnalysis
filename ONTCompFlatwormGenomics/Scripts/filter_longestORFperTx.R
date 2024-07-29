library(tidyverse)

# Simple script to extract the longest ORF per transcripts

args = commandArgs(trailingOnly=TRUE)

gff3_input <- args[1]
gff3_output <- args[2]

cds <- rtracklayer::import(gff3_input)

names(cds) <- paste0("e", str_pad(1:length(cds), width=8, side="left", pad="0"))

cds_uniq <- cds %>% GenomicRanges::mcols() %>%
		as_tibble() %>%
		dplyr::filter(type == "mRNA") %>%
		dplyr::select(Name, ORF_Len) %>%
		mutate(tx = gsub(".p",";",Name, fixed =T)) %>%
		separate(tx, into="tx", sep=";", extra="drop") %>%
		group_by(tx) %>%
		slice_max(ORF_Len, n=1, with_ties = FALSE)

mrna <- cds[cds$type == "mRNA" & cds$Name %in% cds_uniq$Name]
feats <- cds[cds$type != "mRNA" & unlist(cds$Parent) %in% cds_uniq$Name]

cds_out <- c(mrna, feats)
cds_out <- cds_out[order(names(cds_out))]

rtracklayer::export(cds_out, gff3_output)
