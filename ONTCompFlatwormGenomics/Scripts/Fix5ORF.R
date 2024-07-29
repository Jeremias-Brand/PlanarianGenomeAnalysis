library(tidyverse)
library(Biostrings)
library(BSgenome)
library(pbapply)

args = commandArgs(trailingOnly=TRUE)

input_gff3 <- args[1]
genome_file <- args[2]
output_gff3 <- args[3]

genome <- readDNAStringSet(genome_file)

# Import Raw GFF3
gff3 <- rtracklayer::import(input_gff3)

# Freeze row order

names(gff3) <- str_pad(1:length(gff3), width=8, side="left", pad="0")

# Split for processing must be done on ID (cds.h1SMcT0000127.1.p1), on CDS only...

partial5p <- mcols(gff3) %>% 
		as_tibble() %>% 
		filter(ORF_Type == "5prime_partial") %>%
		pull(ID)

partial5p <- paste0("cds.",partial5p)

gff3_2fix <- gff3[gff3$type == "CDS"]
gff3_2fix <- split(gff3_2fix, gff3_2fix$ID)

gff3_2fix <- gff3_2fix[names(gff3_2fix) %in% partial5p]

gff3$NewCDSStart <- NA

# Core program

fix5pORF <- function(i){

	goi <- gff3_2fix[[i]]

	# Consider first exon, according to strandedness
	strand = as.vector(GenomicRanges::strand(goi)[1])

	if (strand == "+") {
		first_exon <- goi[which.min(start(goi))]
	}

	if (strand == "-") {
		first_exon <- goi[which.max(end(goi))]
	}

	if (strand == ".") {
		error("Something wrong with strandedness!")
	}

	# Split CDS exon into codons
	seq <- unlist(strsplit(gsub("(.{3})", "\\1 ", getSeq(genome, first_exon)), " "))

	# Pick first AUG codon, if any
	start_codon <- min(which(seq == "ATG"))
	cds_start_pos = (3*(start_codon-1)+1)

	# Convert codon number to position within CDS
	if (strand == "+") {
		cds_start_pos_genom = GenomicRanges::start(first_exon) + cds_start_pos - 1
		first_exon$NewCDSStart <- cds_start_pos_genom
	}

	if (strand == "-") {
		cds_start_pos_genom = GenomicRanges::end(first_exon) - cds_start_pos + 1
		first_exon$NewCDSStart <- cds_start_pos_genom		
	}

	return(first_exon)
}

# Parallel execution

library(parallel)
cl <- makeCluster(3)
clusterExport(cl, c("fix5pORF", "gff3_2fix", "getSeq", "genome"))

corrected <- pblapply(partial5p, fix5pORF, cl=cl) %>% 
		GRangesList() %>% 
		unlist()

stopCluster(cl)

# Purge gff from uncorrected lines

gff3 <- gff3[! names(gff3) %in% names(corrected)]

# Replace corrected rows

gff3 <- c(gff3, corrected)

gff3 <- gff3[order(names(gff3))]

# Calculate absolute AUG distance from tx 5' end

gff3$DeltaTx = abs( gff3$NewCDSStart - start(gff3) )

# Calculate 5'UTR lenght for transcripts with complete ORFs
# (Ground reference)

library(GenomicFeatures)

complete_tx <- gff3[gff3$type=="mRNA" & gff3$ORF_Type == "complete"]$ID

txdb <- makeTxDbFromGFF(input_gff3)
utrl <- transcriptLengths(txdb, with.utr5_len=T) %>% as_tibble()

utrl <-  utrl %>% filter(tx_name %in% complete_tx) 

# Estimate threshold for maximum 5'utr lenght

library(mixtools)
lutrl <- log2(utrl$utr5_len)
lutrl <- lutrl[is.finite(lutrl)]
gmm <- normalmixEM(lutrl)

thresh_len = 2^(gmm$mu[which.min(gmm$mu)] + 3*gmm$sigma[which.min(gmm$mu)])

pdf(paste0(input_gff3,".utr5p_lenght_distribution.pdf"), width=4, height=2.5)
plot(gmm, which=2, breaks=200)
abline(v=log2(thresh_len), lty=3, col="red")
dev.off()

# Set threshold for considering a true start sitethresh_len

# before_fix <- gff3[!is.na(gff3$ORF_Type)]$ORF_Type %>% table()

# Fix start positions and 

gff3[is.na(gff3$DeltaTx)]$DeltaTx <- Inf

start(gff3[gff3$DeltaTx < thresh_len & strand(gff3) == "+"]) <- gff3[gff3$DeltaTx < thresh_len & strand(gff3) == "+"]$NewCDSStart
end(gff3[gff3$DeltaTx < thresh_len & strand(gff3) == "-"]) <- gff3[gff3$DeltaTx < thresh_len & strand(gff3) == "-"]$NewCDSStart

# Rename ORF_types

ORFtype_to_rename <- gsub("cds.", "", gff3[gff3$DeltaTx < thresh_len]$ID, fixed=T)

gff3[gff3$ID %in% ORFtype_to_rename & gff3$type == "mRNA"]$ORF_Type <- "complete"

# after_fix <- gff3[!is.na(gff3$ORF_Type)]$ORF_Type %>% table()

gff3$Name[gff3$type != "mRNA"] <- NA 
gff3$NewCDSStart <- NULL
gff3$DeltaTx <- NULL

gff3 <- gff3[gff3$type != "three_prime_UTR"]
gff3 <- gff3[gff3$type != "five_prime_UTR"]

# Export

rtracklayer::export(gff3, output_gff3)
write_tsv(as_tibble(ORFtype_to_rename), paste0(output_gff3, ".fixed_ORFs.tsv"))
