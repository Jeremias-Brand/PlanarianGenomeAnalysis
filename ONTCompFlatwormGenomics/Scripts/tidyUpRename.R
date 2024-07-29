library(tidyverse)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
input1 <- args[1]
input2 <- args[2]
output <- args[3]
rootname <- args[4]

# Import and merge positive/negative strands

pos <- rtracklayer::import(input1)
pos$locus[!is.na(pos$locus)] <- paste0(pos$locus[!is.na(pos$locus)],"_pos")

neg <- rtracklayer::import(input2)
neg$locus[!is.na(neg$locus)] <- paste0(neg$locus[!is.na(neg$locus)],"_neg")

annot <- c(pos, neg)
rm(pos, neg)

annot$gene_id[annot$type == "exon"] <- NA

# Order loci by position

annot$orderpt <- BiocGenerics::start(annot)
annot$orderchr <- as.vector(seqnames(annot))
annot$scaffold <- grepl("scaffold", annot$orderchr)

annot$orderchr[annot$scaffold] <- tibble(orderchr=annot$orderchr[annot$scaffold]) %>%
			separate(orderchr, into=c("V1","V2", "V3")) %>%
			mutate(V2=str_pad(V2,3,"left","0")) %>%
			unite(orderchr, V1, V2, V3) %>%
			pull(orderchr)

t_order <- mcols(annot[annot$type=="transcript"]) %>%
		as_tibble() %>%
		dplyr::select(orderchr, orderpt, locus) %>%
		na.exclude() %>% 
		group_by(locus) %>%
		summarise(orderchr = head(orderchr, 1), orderpt = min(orderpt)) %>%
		arrange(orderchr, orderpt) %>%
		dplyr::select(-orderchr, -orderpt)

t_order$new_gene_id <- paste0(rootname, "G", str_pad(1:nrow(t_order), width=7, side="left", pad="0")) 

annot$gene_id <- annot$locus
annot$locus <- NULL
annot$orderchr <- NULL
annot$orderpt <- NULL
annot$scaffold <- NULL

# Rename loci and transcripts

name_t <- mcols(annot) %>%
		as_tibble() %>%
		dplyr::select(transcript_id, gene_id) %>%
		na.exclude() %>%
		unique() %>%
		left_join(t_order, by=c("gene_id"="locus"))
	
name_t <- split(name_t, name_t$new_gene_id)
			
for (n in names(name_t)){
	name_t[[n]]$new_transcript_id <- paste0(n, ".", 1:nrow(name_t[[n]]))
}

name_t <- bind_rows(name_t) %>%
		 mutate(new_transcript_id = gsub(paste0(rootname,"G"), paste0(rootname,"T"), new_transcript_id))
		 
# Apply new names to the annotations

annot$new_gene_id <- name_t$new_gene_id[match(annot$gene_id, name_t$gene_id)]
annot$new_transcript_id <- name_t$new_transcript_id[match(annot$transcript_id, name_t$transcript_id)]

annot$gene_id <- annot$new_gene_id
annot$transcript_id <- annot$new_transcript_id
annot$new_gene_id <- NULL
annot$new_transcript_id <- NULL
annot$lineorder <- NULL

tx2g <- mcols(annot) %>% 
		as_tibble() %>% 
		dplyr::select(transcript_id) %>%
		unique() %>%
		separate(transcript_id, into="gene_id", extra="drop", remove=F) %>%
		mutate(gene_id=gsub(paste0(rootname,"T"), paste0(rootname,"G"),gene_id))
		
annot$gene_id <- tx2g$gene_id[match(annot$transcript_id, tx2g$transcript_id)]
annot$score <- NA

rtracklayer::export(annot, output)
