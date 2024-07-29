library(tidyverse)
library(Biostrings)
library(pbapply)
library(parallel)
library(GenomicFeatures)

# This script purges transcript fragments 
# (encoding portions of the full ORF)

args = commandArgs(trailingOnly=TRUE)

#GTF_INPUT <- "Snov.mix_annotation-V5MF22CN.gtf"
#CDS_PEP <- "Snov.mix_annotation-V5MF22CN.CDS.pep"
#GTF_OUTPUT_coding <- "Snov.mix_annotation-V5MF22CN_FiltIso.gtf"
#GTF_OUTPUT_nc <- "Snov.mix_annotation-V5MF22CN_nc.gtf"
#cdhit_path <- "/home/luca/bin/cd-hit-v4.8.1-2019-0228/"

GTF_INPUT <- args[1]
CDS_PEP <- args[2]
GTF_OUTPUT_coding <- args[3]
GTF_OUTPUT_nc <- args[4]
cdhit_path <- args[5]

# Consider only multi-exonic transcripts

txdb <- makeTxDbFromGFF(GTF_INPUT)
tx <- transcriptLengths(txdb)

nExon <- exonsBy(txdb, by="tx", use.names=T) %>%
		elementNROWS() %>% 
		enframe() %>%
		dplyr::rename(tx = name, nEx = value)
		
multiExonTx <- filter(nExon, nEx>1) %>% pull(tx)

# Import raw ORF content info from TransDecoder

peptides <- readAAStringSet(CDS_PEP)

peptides <- as.vector(peptides) %>% 
		enframe() %>%
		dplyr::rename(seq = value, ID = name) %>%
		mutate(ID = gsub(" ",";", ID, fixed=T)) %>%
		separate(ID, into="ID", extra="drop", sep=";") %>%
		mutate(txtmp = gsub(".",";", ID, fixed=T)) %>%
		separate(txtmp, into=c("Gene","tx"), extra="drop", sep=";") %>%
		mutate(tx = paste0(Gene, ".",tx)) %>%
		mutate(Gene = gsub("T","G", Gene))

# Consider only multiexonic transcripts for pruning

pepByGene <- filter(peptides, tx %in% multiExonTx) %>%
		group_by(Gene) %>%
		group_split(.keep=T)

# es: h1SMG0011538 foo <- pepByGene[[4621]]

# Run genewise ORF clustering with CD-HIT (90% identity)
# and merge cluster information

cdhit <- function(x){ 

	rndl <- paste0(sample(c("0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f"), 12), collapse="") # Unique tag ID for tmp files
	aas <- Biostrings::AAStringSet(x$seq) 
	names(aas) <- x$ID
	rtracklayer::export(aas, paste0("/dev/shm/", rndl,".fasta"))
	foo <- system(paste0(cdhit_path,"/cd-hit -i /dev/shm/",rndl,".fasta -o /dev/shm/",rndl,"_cd.fasta -c 0.9"),intern=TRUE)
	foo <- system(paste0(cdhit_path,"/clstr2txt.pl /dev/shm/",rndl,"_cd.fasta.clstr > /dev/shm/",rndl,"_cd.txt"),intern=TRUE)

	suppressMessages(clstr <- readr::read_tsv(paste0("/dev/shm/", rndl, "_cd.txt"), progress=F) %>% dplyr::select(id, clstr, length, clstr_iden, clstr_cov, clstr_rep))

	system(paste0("rm /dev/shm/", rndl,"*"))

	dplyr::left_join(x, clstr, by=c("ID"="id")) %>%
	dplyr::mutate(clstr_cov = as.numeric(gsub("%","", clstr_cov, fixed=T)),
		clstr_iden = as.numeric(gsub("%","", clstr_iden, fixed=T))) %>%
	return()
}

# Filter out transcripts producing ORFs with a coverage lower than 85%
# (referred to the cluster leader)
# Discard tx only if all p* generated from that transcript 
# have to be discarded, i.e. it has max(clstr_cov)<85, otherwise keep:

dismiss <- function(x){
		cdhit(x) %>% 
			dplyr::mutate(ID = gsub(".p",";", ID, fixed=T)) %>%
			tidyr::separate(ID, into="tx", extra="drop", remove=F, sep=";") %>%
			dplyr::mutate(ID = gsub(";",".p", ID, fixed=T)) %>%
			dplyr::select(-seq, -Gene) %>%	
			dplyr::group_by(tx) %>%
			dplyr::summarise(max_clstr_cov = max(clstr_cov)) %>%
			dplyr::filter(max_clstr_cov < 85) %>%
			dplyr::pull(tx) %>%
			return()
		}

cl <- makeCluster(18)
clusterExport(cl, c("cdhit", "cdhit_path", "dismiss", "%>%"))

pbo = pboptions(type="txt")

tx_discarded <- pblapply(pepByGene, dismiss, cl=cl) %>% unlist()

stopCluster(cl)

annot <- rtracklayer::import(GTF_INPUT)

# Filter out discarded transcripts

filter_annot <- annot[! annot$transcript_id %in% tx_discarded]

# Extract non-coding transcript into a putative ncRNA collection

nc_transcripts <- filter_annot[! filter_annot$transcript_id %in% peptides$tx]
coding_transcripts <- filter_annot[filter_annot$transcript_id %in% peptides$tx]

# Export final files

rtracklayer::export(coding_transcripts, GTF_OUTPUT_coding)
rtracklayer::export(nc_transcripts, GTF_OUTPUT_nc)

