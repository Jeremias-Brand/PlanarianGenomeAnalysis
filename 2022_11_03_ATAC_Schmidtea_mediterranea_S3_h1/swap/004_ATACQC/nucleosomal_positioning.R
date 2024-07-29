#!/usr/bin/env Rscript

library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(ATACseqQC)
library(ChIPpeakAnno)
library(VplotR)
library(BiocManager)
library(BSgenome.Smediterranea.MPINAT.S3h1)

load(file = "ATACseqQC.RData")
path_to_split <- "split_bam"

args = commandArgs(trailingOnly=TRUE)
string_pre1 = "NFR"
string_pre2 = "mononucleosomal"
string_pre3 = "dinucleosomal"
string_pre4 = "trinucleosomal"
string_suff= "sorted.bam"
 
# concatenate two strings using separator
NFR = paste(string_pre1, args[1], string_suff, sep = "_")
mononucleosomal = paste(string_pre2, args[1], string_suff, sep = "_")
dinucleosomal = paste(string_pre3, args[1], string_suff, sep = "_")
trinucleosomal = paste(string_pre4, args[1], string_suff, sep = "_")


##read in bam-files
split_bamfiles <- file.path(path_to_split,
                            c(NFR,
                              mononucleosomal,
                              dinucleosomal,
                              trinucleosomal
				))

split_bamfiles.labels <- gsub(".bam", "", basename(split_bamfiles))

## estimate the library size for normalization
librarySize <- estLibSize(split_bamfiles)

NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(split_bamfiles,
                          TSS=TSS_smed,
                          librarySize=librarySize,
                          seqlev=seqinformation_txdb@seqnames,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)

out_Nucl <- featureAlignedDistribution(sigs, reCenterPeaks(TSS_smed, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")

out_Nucl <- apply(out_Nucl, 2, range01)

pdf(file = args[2], width = 12,height = 9.5)

matplot(out_Nucl, type="l", xaxt="n", 
        main = paste0("Nucleosomal-free and mononucleosomal signal at TSS for sample ", split_bamfiles.labels),
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)

abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
legend("topright", colnames(out_Nucl),col=seq_len(ncol(out_Nucl)),cex=.8,fill=seq_len(ncol(out_Nucl)))
dev.off()
