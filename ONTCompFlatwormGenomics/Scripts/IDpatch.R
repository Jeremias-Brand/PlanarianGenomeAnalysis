library(tidyverse)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

input_gff3 <- args[1]
output_gff3 <- args[2]

# Fix small GFF3 format inconsistencies
# (add gene names, rename source, fix fields for symmetry)

gr <- rtracklayer::import(input_gff3)

names(gr) <- paste0("e", 1:length(gr))

edit_rows <- gr[gr$type %in% c("transcript", "mRNA")]

gr <-  gr[! gr$type %in% c("transcript", "mRNA")]

edit_rows$Name <- edit_rows$ID

# Add Gene ID
edit_rows$geneID <- edit_rows$ID %>% 
			as_tibble() %>%
			separate(value, into="value", extra="drop") %>%
			mutate(value = gsub("T","G", value)) %>%
			pull(value)

edit_rows$Parent <- edit_rows$geneID
edit_rows$type <- "transcript"

genes <- as_tibble(edit_rows) %>%
		group_by(geneID) %>%
		summarise(seqnames=dplyr::first(seqnames), start=min(start), end=max(end), strand=dplyr::first(strand))
	
genes <- GRanges(seqnames=genes$seqnames, IRanges(genes$start, genes$end), strand=genes$strand, type="gene", ID=genes$geneID, Name=genes$geneID)

gr <- c(gr, edit_rows, genes)
gr$source <- "ONThybrid"

rtracklayer::export(gr, con=output_gff3)
