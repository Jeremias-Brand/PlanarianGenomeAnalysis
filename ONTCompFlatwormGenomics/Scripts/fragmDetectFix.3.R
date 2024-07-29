library(tidyverse)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)

# Legend
# M : mixed, before FLAIR
# f : after FLAIR

M_GTF = args[1]
M_CDS_GFF3 = args[2]
M_NAMES = args[3]

f_GTF = args[4]
f_CDS_GFF3 = args[5]
f_NAMES = args[6]

TX_NAME_MAP = args[7]
OUT_ANN = args[8]

# #####################################################################################
# Input ORFs
# #####################################################################################

ORF_M <- rtracklayer::import(M_CDS_GFF3) %>% 
	mcols() %>% 
	as_tibble() %>%
	filter(type == "mRNA") %>%
	mutate(Parent=unlist(Parent),
		ORF_Len = as.numeric(ORF_Len))

ORF_M_genes <- read_delim(M_NAMES, col_names=c("tx","gene"), delim=" ") %>%
		mutate(tx = gsub(">","",tx, fixed=T))

ORF_M  <- ORF_M %>% 
	dplyr::rename(len = ORF_Len) %>%
	dplyr::select(ID, len) %>%
	mutate(ID = gsub(".p",";", ID, fixed=T)) %>%
	separate(ID, into=c("ID"), sep=";", extra="drop") %>%
	left_join(ORF_M_genes, by=c("ID"="tx"))

#> ORF_M
## A tibble: 89,927 × 3
#   ID          len gene      
#   <chr>     <dbl> <chr>     
# 1 H.1.1        78 H.1       
# 2 H.10.1       76 MSTRG.6   
# 3 H.100.1      99 MSTRG.56  
# 4 H.1000.1     90 H.1000    
# 5 H.1000.2     72 H.1000    

# #####################################################################################
# FLAIR ORFs
# #####################################################################################

ORF_F <- rtracklayer::import(f_CDS_GFF3) %>% 
	mcols() %>% 
	as_tibble() %>%
	filter(type == "mRNA") %>%
	mutate(Parent=unlist(Parent),
	ORF_Len = as.numeric(ORF_Len))
	
# Remap "chr:XXXX or XXXX-XXXX genes" to V5M identifiers 

transcript_map <- read_tsv(TX_NAME_MAP) %>%
			filter(class_code %in% c("=","c")) %>%
			dplyr::select(qry_id, ref_gene_id)

#> transcript_map
## A tibble: 77,042 × 2
#   qry_id                                        ref_gene_id
#   <chr>                                         <chr>      
# 1 122:1698|24d5465b-f9a9-4f71-bde6-894cf0f25cb6 MSTRG.5    
# 2 130:1973|fa89de4f-585b-4b6c-ac81-76116c094724 MSTRG.8    
# 3 131:927|caca0d88-d4a7-4352-aa53-d261b6dc0c7f  MSTRG.1    
# 4 249093f4-d339-42dd-b8a6-1ce36466c8fd          MSTRG.15   
# 5 c892f16e-391f-4163-99fd-e368d4c2ae85          MSTRG.3    

ORF_F_genes_raw <- read_delim(f_NAMES, col_names=c("tx","gene"), delim=" ") %>%
		mutate(tx = gsub(">","",tx, fixed=T))

ORF_F_genes_1 <- ORF_F_genes_raw %>% 
		filter(grepl("MSTRG|H.",gene)) %>%
		mutate(XX = F)

#> ORF_F_genes_1
## A tibble: 98,639 × 3
#   tx                                            gene        XX   
#   <chr>                                         <chr>       <lgl>
# 1 H.16189.2                                     MSTRG.10520 FALSE
# 2 H.16189.5                                     MSTRG.10520 FALSE
# 3 H.16189.3                                     MSTRG.10520 FALSE
# 4 H.14534.2                                     MSTRG.9502  FALSE
# 5 140:2430|992d982d-64b7-4123-9139-ecc1d717751f MSTRG.9502  FALSE
# 6 H.14534.3                                     MSTRG.9502  FALSE
# 7 123:1241|64f5254e-61bf-4f6c-a95e-6310b9bd2ee2 MSTRG.9502  FALSE

ORF_F_genes_2 <- ORF_F_genes_raw %>% 
		filter(!grepl("MSTRG|H.",gene)) %>% 
		left_join(transcript_map, by=c("tx"="qry_id")) %>% 
		na.exclude() %>% # => recover only chr:XXXX or XX:XXXX genes overlapping V5M annotation
		dplyr::select(-gene) %>%
		dplyr::rename(gene = ref_gene_id) %>%
		mutate(XX = T)

ORF_F_genes <- bind_rows(ORF_F_genes_1, ORF_F_genes_2)
rm(ORF_F_genes_1, ORF_F_genes_2)

ORF_F  <- ORF_F %>% 
	dplyr::rename(len = ORF_Len) %>%
	dplyr::select(ID, len) %>%
	mutate(ID = gsub(".p",";", ID, fixed=T)) %>%
	separate(ID, into=c("ID"), sep=";", extra="drop") %>%
	left_join(ORF_F_genes, by=c("ID"="tx")) %>%
	na.exclude()

# #####################################################################################
# Compare side by side (sbs)
# #####################################################################################

ORF_M_max <- dplyr::select(ORF_M, -ID) %>%
		group_by(gene) %>%
		summarise(len=max(len))

ORF_F_max <- dplyr::select(ORF_F, -ID) %>%
		group_by(gene) %>%
		summarise(len=max(len))

sbs <- full_join(ORF_M_max, ORF_F_max, by="gene")
sbs[is.na(sbs)]<-0

#> sbs
## A tibble: 37,813 × 3
#   gene    len.x len.y
#   <chr>   <dbl> <dbl>
# 1 H.1        78     0
# 2 H.1000     90     0
# 3 H.10001    69     0
# 4 H.10006   350   350
# 5 H.10010   162     0
# 6 H.10011    82    72
# 7 H.10059   125     0

## Find which transcripts have been fragmented after FLAIR refinement

recov_genes <- filter(sbs, len.x-len.y>20, len.x>120) %>%  # Chosen in order to maximise BUSCOs captured
			pull(gene)
			
recov_tx <- group_by(ORF_M, gene) %>% 
		slice_max(len, n=1, with_ties=F) %>%
		filter(gene %in% recov_genes) %>%
		pull(ID)

# #####################################################################################
# Import and process tx annotations
# #####################################################################################

# Import Source annotation

M <- rtracklayer::import(M_GTF)

# Import Destination (FLAIR) annotation

fl <- rtracklayer::import(f_GTF)

# Remove all "chr:XXXX or XXXX-XXXX" tx

fl0 <- fl[grepl("MSTRG|H.", fl$gene_id)]

# Insert back chr:XXXX or XXXX-XXXX overlapping with V5M genes
# That is: keep all tx that have been reassigned by FLAIR but had a similar tx in V5M
# (i.e. those chr:XXXX in ORF_F_genes)

f_recover <- fl[fl$transcript_id %in% filter(ORF_F_genes, XX)$tx]
f_recover$transcript_id <- paste0(f_recover$transcript_id, "_recf") # Append suffix to tx_id

# Insert back transcripts with longest ORF from V5M

M_recover <- M[M$transcript_id %in% recov_tx]
M_recover$transcript_id <- paste0(M_recover$transcript_id, "_recM") # Append suffix to tx_id

fl2 <- c(fl0, M_recover, f_recover)

# #####################################################################################
# Sort and Export 
# #####################################################################################

fl2 <- sortSeqlevels(fl2)
fl2 <- sort(fl2)

fl2$score <- NA

rtracklayer::export(fl2, OUT_ANN)

