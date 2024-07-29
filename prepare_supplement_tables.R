################################################################################

## collection and preparation of the output tables

################################################################################

library(writexl)
library(readxl)
library(tidyverse)
library(GenomicRanges)
library(data.table)
# Narrow peak header
# We slightly modify the standard naming because we use peak to refere to the entier called region and summit to what is called peak in the narrowPeak format.
narrow_header <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'summit')
################################################################################
## Peaks that were condidates for conservation with Schmidtea mediterranea

################################################################################
## ChIP peak files
## these represent the entier signal of the unified peak sets in each species
## all of them were used in the assessment of conservation because in a few the value for the summit is 0. 
chip_files = list(
  wtH3K4me3 = './final_putative_enhancers/chr_wt_H3K4me3_peaks.narrowPeak',
  wtH3K27ac = './final_putative_enhancers/chr_wt_H3K27ac_peaks.narrowPeak'
)

chip <- lapply(chip_files,
function(x) {
df <- read.table(x)
header <- c( "chrom", "chromStart", "chromEnd",
  "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
names(df) <- header
    df$peak = paste(df$chrom, df$chromStart, df$chromEnd, sep = '_')
    df <- dplyr::select(df, peak, everything())
    return(df)
  }
)

# MACS2 narrow peak files
narrow_files = list(
  schMedS3h1 = './final_putative_enhancers/chr_Smed_combined_peaks.narrowPeak',
  schPol2 = './final_putative_enhancers/chr_Spol_combined_peaks.narrowPeak',
  schLug1 = './final_putative_enhancers/chr_Slug_combined_peaks.narrowPeak',
  schNov1 = './final_putative_enhancers/chr_Snov_combined_peaks.narrowPeak'
)

# CLOSEST GENE ANNOTATION
closest_files = list(
  schMedS3h1 = './peak_annotation/schMedS3h1_by_class_closest.tsv',
  schPol2 = './peak_annotation/schPol2_closest.tsv',
  schLug1 = './peak_annotation/schLug1_closest.tsv',
  schNov1 = './peak_annotation/schNov1_closest.tsv'
)

# species names for interation and checking
sps <- c(
  'schMedS3h1',
  'schPol2',
  'schLug1',
  'schNov1'
)

# checking that the naming is correct
# Function to check if the order of names in a list matches the sps vector
check_order <- function(lst, sps) {
  all(names(lst) == sps)
}
# abort if not in the same order
stopifnot(exprs = {
  check_order(narrow_files, sps)
  check_order(closest_files, sps)
})


narrow <- lapply(narrow_files, function(x) {
  df <- read.table(x)
  names(df) = narrow_header
  return(df)
})

closest <- lapply(closest_files, function(x) {
  df <- read.table(x, header = TRUE) |>
  dplyr::select(-strand) |>
  # narrow peaks is zero indexed 
  mutate(start = start - 1)
  return(df)
})

merged <- purrr::map2(closest, narrow, full_join,  by = c(
  'seqnames' = 'chrom', 
  'start' = 'chromStart',
  'end' = 'chromEnd'
  ))

################################################################################

## conservation data

################################################################################
cons <- lapply(list('conserved_elements/smed_base_full_conservation_schPol2_coordinates.tsv',
            'conserved_elements/smed_base_full_conservation_schLug1_coordinates.tsv',
            'conserved_elements/smed_base_full_conservation_schNov1_coordinates.tsv'),
       read.table, header = TRUE)
names(cons) <- c(
  'schPol2_WGA_schMedS3h1',
  'schLug1_WGA_schMedS3h1',
  'schNov1_WGA_schMedS3h1'
  )

################################################################################

## intersecting the WGA results with the peaks in each species (again)

################################################################################
# Check if the length of 'cons' and 'merged' is sufficient
if(length(cons) < 3 || length(merged) < 4) {
  stop("Error: 'cons' or 'merged' does not have the expected number of elements.")
}

for (k in 1:3) {
  print(paste('Iteration', k))
  liftover = cons[[k]]
  atac = merged[[k + 1]]

  # Check if required columns exist
  if(!("peak" %in% names(atac)) || !("smed_peak" %in% names(liftover))) {
    stop("Error: Required columns not found in data frames.")
  }

  # Convert the first data frame to a GRanges object
  # QUERY
  gr1 <- makeGRangesFromDataFrame(
    liftover,
    keep.extra.columns=TRUE, seqnames.field="seqnames", start.field="start", end.field="end", strand = 'strand' 
  )
  # Convert the second data frame to a GRanges object
  # SUBJECT
  gr2 <- makeGRangesFromDataFrame(
    atac,
    keep.extra.columns=TRUE, seqnames.field="seqnames", start.field="start", end.field="end", strand = 'strand' 
  )

  ################################################################################
  # Find overlap of each liftover with the ATAC-seq peaks in that species.
print(paste('Find overlap of each liftover with the ATAC-seq peaks in that species.'))
  overlaps <- findOverlaps(gr1, gr2)
  # add a results column
  liftover$overlapping_atac_peak <- NA

  # annotate liftover based on these overlaps creating a new row for each overlap
  for (i in unique(queryHits(overlaps))){
    # get the matching entries
    hits_id = subjectHits(overlaps)[queryHits(overlaps) == i] 
    # annotate the query data frame with the subject hit results 
    liftover[i,'overlapping_atac_peak'] = paste0(unique(atac[hits_id,'peak']), collapse = ';')
  } 

  ################################################################################
  # Find overlap of each ATAC-seq peaks with a Liftover peak.
  print(paste('Find overlap of each ATAC-seq peaks with a Liftover peak.'))
  overlaps <- findOverlaps(gr2, gr1)
  # add a results column
  atac$smed_peak_liftover <- NA

  # annotate liftover based on these overlaps creating a new row for each overlap
  for (i in unique(queryHits(overlaps))){
    # get the matching entries
    hits_id = subjectHits(overlaps)[queryHits(overlaps) == i] 
    # annotate the query data frame with the subject hit results 
    atac[i,'smed_peak_liftover'] = paste0(unique(liftover[hits_id,'smed_peak']), collapse = ';')
  } 
  # assign the transformed data.frames back to the originals
  print(paste('Done, assigning output.'))
  cons[[k]] = liftover
  merged[[k + 1]] = atac
}

################################################################################

## Adding orthology data

gs <- fread('./genespace/results/combBed.txt')
# filter the non-duplicated high confidence conserved peaks
for (NAME in names(merged)) {
  # Extract dataframe from 'merged' list and create a key by removing numerical suffix from 'transcriptId'
  df = merged[[NAME]] |>
    mutate(key = gsub('\\.[0-9]+$', '', transcriptId))

  # Filter the 'gs' dataframe for the current genome, create a key, and select relevant columns
  og_tbl <- gs |>
    filter(genome == NAME) |>
    dplyr::select(id, og) |>
    mutate(key = gsub('\\.[0-9]+$', '', id)) |>
    dplyr::select(key, og)

  # Perform a left join of 'df' with 'og_tbl' using the created key
  res = left_join(df, og_tbl, by = 'key')

  # Assign the result back to the 'merged' list
  merged[[NAME]] = res
}

# smed annotated peaks
NAME = 'schPol2'
s1 <- merged[[1]]

for (k in 1:3) {
NAME = names(merged)[k + 1]
x <- cons[[k]] |> 
filter(!is.na(overlapping_atac_peak) & summit_OL_count == 1) |>
  dplyr::select(overlapping_atac_peak, smed_peak) |>
  unique()

# there are some cases where a S. mediterranea peak liftover overlaps with two different peaks
# we exclude these many-to-many relationships
name_duplicates <- x |> group_by(overlapping_atac_peak) |> summarise(N = n()) |> filter(N > 1) |> pull(overlapping_atac_peak)
smed_duplicates <- x |> group_by(smed_peak) |> summarise(N = n()) |> filter(N > 1) |> pull(smed_peak)

x <- x |> filter(!overlapping_atac_peak %in% name_duplicates & ! smed_peak %in% smed_duplicates)

# add the closest gene information
xx <- left_join(x, dplyr::select(merged[[k + 1]], peak, geneId, seqnames, geneStart, geneEnd, og), by = c('overlapping_atac_peak' = 'peak'))
xx <- xx |>
dplyr::rename_with(~ paste0(NAME, "_", .), names(xx)[names(xx) %in% names(s1)]) |>
dplyr::rename_with(~ paste0(NAME, "_peak"), overlapping_atac_peak)

s1 <- left_join(s1, xx, by = c('peak' = 'smed_peak'))
}

s1 <- s1 |>
mutate(same_og = ifelse(og == schPol2_og & og == schNov1_og & og == schLug1_og, 'YES', 'NO'))

#################################################

# Same filtering but less duplicate removal.
## merger including all duplications and 
df <- merged[[1]]
linked_smed = list()
# loop through each species
for (k in 1:3) {
NAME = names(merged)[k + 1]
con_df <- cons[[k]] |> 
filter(!is.na(overlapping_atac_peak) & summit_OL_count == 1) |>
  dplyr::select(overlapping_atac_peak, smed_peak) |>
  unique()

# add the closest gene information
con_df2 <- left_join(con_df, dplyr::select(merged[[k + 1]], peak, geneId, seqnames, geneStart, geneEnd, og), by = c('overlapping_atac_peak' = 'peak'))
con_df2 <- con_df2 |>
dplyr::rename_with(~ paste0(NAME, "_", .), names(con_df2)[names(con_df2) %in% names(df)]) |>
dplyr::rename_with(~ paste0(NAME, "_peak"), overlapping_atac_peak)

linked_smed[[NAME]] <- left_join(df, con_df2, by = c('peak' = 'smed_peak'))
}

write_xlsx(linked_smed, path = './supporting_information/conserved_peaks_per_species.xlsx')


###################################

## read block coords

blk <- fread('./genespace/results/blkCoords.txt')
###################################

## OUTPUT TABLE FOMATTING

####################################

# CREATE SUPPORTING FILE

SM_list <- list(
  TabS1_schMedS3h1_ATAC = s1,
  TabS2_schMedS3h1_wtH3K4me3 = chip[[1]],
  TabS3_schMedS3h1_wtH3K27ac = chip[[2]],
  TabS7_schPol2_ATAC = merged$schPol2, 
  TabS8_schNov1_ATAC = merged$schNov1, 
  TabS9_schLug1_ATAC = merged$schLug1, 
  TabS10_schPol2_WGA = cons$schPol2_WGA_schMedS3h1, 
  TabS11_schNov1_WGA = cons$schNov1_WGA_schMedS3h1, 
  TabS12_schLug1_WGA = cons$schLug1_WGA_schMedS3h1, 
  TabS13_GENESPACE = gs,
  TabS14_synteny_blocks = blk
)

####################################

# Write table

write_xlsx(SM_list, './supporting_information/Additional\ File\ 2.xlsx')


###

## Summary tables for Additional File 1

##

df <- readxl::read_excel('./supporting_information/Additional\ File\ 2.xlsx', sheet = 1)

df |> 
group_by(conservation_class, element_type)  |> 
mutate(conservation_class = factor(conservation_class, levels = c(
  'highly_conserved', 'partially_conserved', 'not_conserved'
)),
element_type = factor(element_type, levels = c(
  'putative_promoter', 'putative_enhancer', 'uncharacterized'
)) ) |> 
summarise(N = n()) |> 
pivot_wider(names_from = element_type, values_from = N)  |> 
write_xlsx('./supporting_information/element_type_conservation_class_N.xlsx')


df |> 
group_by(simple_annotation, element_type)  |> 
mutate(simple_annotation = factor(simple_annotation, levels = c(
  '<1kb TSS', '1-2kb TSS', '2-3kb TSS',
  '5UTR', 'Exon', 'Intron', 'Downstream <=300bp', 'Distal intergenic'
)),
element_type = factor(element_type, levels = c(
  'putative_promoter', 'putative_enhancer', 'uncharacterized'
)) ) |> 
summarise(N = n()) |> 
pivot_wider(names_from = element_type, values_from = N)  |> 
write_xlsx('./supporting_information/element_type_simple_annotation_N.xlsx')


