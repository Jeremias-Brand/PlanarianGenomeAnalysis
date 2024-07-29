library(regioneR)
library(GenomicRanges)
library(readxl)
library(tidyverse)
library(writexl)
################################################################################
permutations = 10000
## we need a BSgenome package for regioneR or just the info from it...

# devtools::install('./SchmedS3_h1_BSgenome/BSgenome.Smediterranea.MPINAT.S3h1')

library(BSgenome.Smediterranea.MPINAT.S3h1)

# Extract chromosome information
chr_info <- seqlengths(BSgenome.Smediterranea.MPINAT.S3h1)

# Create a GRanges object representing each chromosome
genome_gr <- GRanges(seqnames = names(chr_info),
                     ranges = IRanges(start = 1, end = chr_info))
################################################################################
## https://bioconductor.org/packages/release/bioc/vignettes/regioneR/inst/doc/regioneR.html

## Load peak files
file  <- './supporting_information/Additional File 2.xlsx'
atac  <- read_xlsx(file, sheet = "TabS1_schMedS3h1_ATAC")
me  <- read_xlsx(file, sheet = "TabS2_schMedS3h1_wtH3K4me3")
ac  <- read_xlsx(file, sheet = "TabS3_schMedS3h1_wtH3K27ac")


## make genomic ranges
atac_gr <- makeGRangesFromDataFrame(atac, keep.extra.columns = FALSE)
me_gr <- makeGRangesFromDataFrame(me, keep.extra.columns = FALSE)
ac_gr <- makeGRangesFromDataFrame(ac, keep.extra.columns = FALSE)



## seqlevels need to match
seqlevels(atac_gr) <- seqlevels(BSgenome.Smediterranea.MPINAT.S3h1)
seqlevels(me_gr) <- seqlevels(BSgenome.Smediterranea.MPINAT.S3h1)
seqlevels(ac_gr) <- seqlevels(BSgenome.Smediterranea.MPINAT.S3h1)

pt_atac_me <- overlapPermTest(A = atac_gr, B = me_gr, 
 ntimes = permutations,
 force.parallel = TRUE, genome = genome_gr,
                              verbose = TRUE)

pt_atac_me$numOverlaps
# Perform overlap permutation test for ATAC and H3K27ac
pt_atac_ac <- overlapPermTest(A = atac_gr, B = ac_gr, 
                              ntimes = permutations,
                              force.parallel = TRUE, genome = genome_gr,
                              verbose = TRUE)

# Function to extract relevant data from the permTestResultsList object
extractPermTestResults <- function(permTestResult) {
    subResult <- permTestResult$numOverlaps
    data.frame(
        Observed = subResult$observed,
        Mean_Simulated = mean(subResult$permuted),
        ZScore = subResult$observed,
        PValue = subResult$pval,
        NIter = length(subResult$permuted)
    )
}

# Extract results for ATAC vs ME and ATAC vs AC
df_atac_me <- extractPermTestResults(pt_atac_me)
df_atac_ac <- extractPermTestResults(pt_atac_ac)

# Save results to Excel
write_xlsx(list(atac_me = df_atac_me, atac_ac = df_atac_ac), 
           path = "./final_putative_enhancers/schMedS3h1_ChIP_ATAC_overlap_permutation_test.xlsx")

pdf('./final_putative_enhancers/schMedS3h1_ChIP_ATAC_overlap_permutation_test.pdf')
plot(pt_atac_me)
plot(pt_atac_ac)
dev.off()

