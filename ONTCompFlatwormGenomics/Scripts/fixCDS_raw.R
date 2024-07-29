library(tidyverse)

# TransDecoder has a bug in the way it returns the appropriate CDS length and type within gff3 files.
# Here we recover the right information from the fasta headers contained in ".cds" 

args = commandArgs(trailingOnly=TRUE)
input1 <- args[1]	# GFF3 File
input2 <- args[2]	# .cds File
output <- gsub(".rawCDS.gff3", ".CDS.gff3", input1)

pre <- rtracklayer::import(input1)

pre <- pre[pre$type != "gene"]

pre$Name[!is.na(pre$Name)] <- pre$ID[!is.na(pre$Name)]

right_pep <- Biostrings::readDNAStringSet(input2) %>% 
		names() %>% 
		enframe() %>%
		mutate(value = gsub(" ", ";", value, fixed=T)) %>%
		separate(value, into=c("V1","V2","V3"), sep=";", extra="drop") %>%
		dplyr::select(-name) %>%
		mutate(	V2 = gsub("type:","", V2),
			V3 = gsub("len:","", V3) ) %>%
		dplyr::rename(ORF_Type = V2,
				ORF_Len = V3)

pre$ORF_Type <- NA
pre$ORF_Len <- NA

pre$ORF_Type[!is.na(pre$Name)] <- right_pep$ORF_Type[match(pre$Name[!is.na(pre$Name)], right_pep$V1)]
pre$ORF_Len[!is.na(pre$Name)] <- right_pep$ORF_Len[match(pre$Name[!is.na(pre$Name)], right_pep$V1)]

rtracklayer::export(pre, con=output)
