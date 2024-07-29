# This script finds different kinds of Chimaera Events:
#
# 1. Chimaeric Loci 
# (Bundle of physically different but somewhat overlapping transcripts 
# that are assigned to the same Gene)
# This is the most common event (e.g. SMNXG0000203)
#
# 2. Chimaeric Transcripts
# (Transcripts that encode two distinct ORFs 
# but are assembled in a single fused tx )
# Rare but harder to fix (e.g. SMNXG0011702)
#
# 3. Mix of both:(e.g. SMNXG0002453)

library(GenomicFeatures)
library(tidyverse)
library(pbapply)

args = commandArgs(trailingOnly=TRUE)

GTF_INPUT <- args[1]     #mix_annotation-V5MF22CN_FiltIsoN.h1.gtf
GFF_CDS_INPUT <- args[2] #mix_annotation-V5MF22CN_FiltIsoN.h1.CDS.gff3.gz
RAW_CDS_INPUT <- args[3] #mix_annotation-V5MF22CN_FiltIsoN.h1.rawCDS.gff3
GTF_OUTPUT <- args[4]    #mix_annotation-V5MF22CN_FiltIsoN_chfx.h1.gtf

annot <- rtracklayer::import(GTF_INPUT)
annot$Info <- annot$Name

annot$Name[!is.na(annot$Name)] <- annot$ID[!is.na(annot$Name)]

txdb <- makeTxDbFromGRanges(annot)

nExon <- exonsBy(txdb, by="tx", use.names=T) %>%
		elementNROWS() %>% 
		enframe() %>%
		dplyr::rename(tx = name, nEx = value)
		
multiExonTx <- filter(nExon, nEx>1) %>% pull(tx)

# 0. Pre-parse ORF database by length

cds_scores <- mcols(rtracklayer::import(GFF_CDS_INPUT)) %>% 
		as_tibble() %>%
		dplyr::filter(type =="mRNA") %>%
		mutate(ORF_Len = as.numeric(ORF_Len)) %>%
		dplyr::select(Name, ORF_Type, ORF_Len)

good_cds <- filter(cds_scores, ORF_Len > 100) %>% # Should be set by default
		mutate(tx = gsub(".p",";", Name)) %>%
		separate(tx, into="tx", extra="drop", sep=";") %>%
		filter(tx %in% multiExonTx) %>% # Only multi-exonic
		pull(Name)

# 1. by orf (.p*) , take max and min coord (range block)

cds <- rtracklayer::import(RAW_CDS_INPUT)
cds <- cds[! cds$type %in% c("gene")]
cds$Name[cds$type == "mRNA"] <- cds$ID[cds$type == "mRNA"]

cds_rb <- cdsBy(makeTxDbFromGRanges(cds), by="tx", use.names=T) %>% range()

cds_rb <- cds_rb[names(cds_rb) %in% good_cds]

# 2. reduce all range blocks by locus

locus_rb <- unlist(cds_rb)

locus_rb$Name <- tibble(Name=names(locus_rb)) %>%
			separate(Name, into =c("V1","V2","V3")) %>%
			dplyr::rename(Name=V1) %>%
			pull(Name)

locus_rb <- GRangesList(split(locus_rb, locus_rb$Name)) %>% 
		GenomicRanges::reduce() %>%
		unlist()
		
# 3. keep only the loci with > 1 range blocks

chimeras <- tibble(Name = names(locus_rb)) %>% 
		group_by(Name) %>% 
		summarise(n=n()) %>%
		filter(n >1) %>%
		pull(Name)

# length(chimeras)
# [1] 472 Before
# [1] 383 After

enframe(chimeras) %>% dplyr::select(value) %>% write_tsv(paste0(GTF_INPUT,".chimeras.tsv"), col_names=F)

# Curate Chimeric locus topology
# For each locus, map to each transcript its overlapping rangeblock(s)

chim_locus_rb <- locus_rb[names(locus_rb) %in% chimeras]
chim_locus_rb <- split(chim_locus_rb, names(chim_locus_rb))

txlist <- txlist <- exonsBy(txdb, by="tx", use.names=T) %>% unlist()
txlist$Locus <- enframe(names(txlist)) %>% 
			separate(value, into="locus", extra="drop") %>% 
			pull(locus)
txlist$tx <- names(txlist)
txlist <- split(txlist, txlist$Locus)

matchRbLoci <- function(loc){
	t_txlist <- txlist[[loc]]
	t_txlist <- split(t_txlist, t_txlist$tx)

	ola <- findOverlaps(t_txlist, chim_locus_rb[[loc]])

	rb_map <- tibble(tx = names(t_txlist)[ola@from], rb = paste0("rb",ola@to)) %>%
			group_by(tx) %>%
			summarise(rb = first(rb), n=n())

	# If we manage to represent all rangeblocks using non-chimeric transcripts
	# Remove the chimeric transcript(s) as redundant,
	# Otherwise keep it/them and flag it/them

	n_rb <- filter(rb_map, n == 1) %>% pull(rb) %>% unique() %>% length()

	if( n_rb == length(chim_locus_rb[[loc]]) ) {
		tag <- "remove"
	} else {
		tag <- "chim"
	}

	# Tag "Hard" chimeric transcripts
	# A "Hard" chimeric transcript overlaps more than one rangeblock
	rb_map$rb[rb_map$n>1] <- tag
	
	rb_map %>% separate(tx, into = "gene", extra="drop", remove=F) %>%
			mutate(gene = gsub("T","G", gene)) %>%
			unite(geneTo, gene, rb, sep="-", remove=F) %>%
			dplyr::select(-n, -rb) %>%
			return()
}

# > matchRbLoci("h1SMcT0005192")
# A tibble: 5 × 3
#  tx              geneTo               gene         
#  <chr>           <chr>                <chr>        
#1 h1SMcT0005192.1 h1SMcG0005192-rb1    h1SMcG0005192
#2 h1SMcT0005192.2 h1SMcG0005192-rb1    h1SMcG0005192
#3 h1SMcT0005192.3 h1SMcG0005192-rb1    h1SMcG0005192
#4 h1SMcT0005192.4 h1SMcG0005192-remove h1SMcG0005192
#5 h1SMcT0005192.5 h1SMcG0005192-rb2    h1SMcG0005192

chim_locus_map <- pblapply(names(chim_locus_rb), matchRbLoci) %>% 
			bind_rows() 

# Purge dispensable chimaeric transcripts

tx_to_remove <- filter(chim_locus_map, grepl("remove", geneTo)) %>% pull(tx)

chim_locus_map <- filter(chim_locus_map, ! grepl("remove", geneTo)) %>% 
			dplyr::select(tx, geneTo)

annot <- annot[! annot$transcript_id %in% tx_to_remove]

# Rename genes corresponding to fixed/flagged chimaeras

names(annot) <- str_pad(1:length(annot), width=7, side="left", pad="0") #add number to recover order afterwards

annot_noedit <- annot[!annot$transcript_id %in% chim_locus_map$tx]
annot_edit <- annot[annot$transcript_id %in% chim_locus_map$tx]

annot_edit$gene_id <- chim_locus_map$geneTo[match(annot_edit$transcript_id, chim_locus_map$tx)]

annot <- c(annot_edit, annot_noedit)
annot <- annot[order(names(annot))]

# Rename all loci sequentially

rootname = args[5]

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
		dplyr::select(orderchr, orderpt, gene_id) %>%
		na.exclude() %>% 
		group_by(gene_id) %>%
		summarise(orderchr = head(orderchr, 1), orderpt = min(orderpt)) %>%
		arrange(orderchr, orderpt) %>%
		dplyr::select(-orderchr, -orderpt)

t_order$new_gene_id <- paste0(rootname, "G", str_pad(1:nrow(t_order), width=7, side="left", pad="0")) 

annot$orderchr <- NULL
annot$orderpt <- NULL
annot$scaffold <- NULL

# Rename loci and transcripts

name_t <- mcols(annot) %>%
		as_tibble() %>%
		dplyr::select(transcript_id, gene_id) %>%
		na.exclude() %>%
		unique() %>%
		left_join(t_order, by=c("gene_id"))
	
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
names(annot) <- NULL

# Export annotation

rtracklayer::export(annot, con=GTF_OUTPUT)
