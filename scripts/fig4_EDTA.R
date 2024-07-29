# look at genespace output
library(data.table)
library(tidyverse)
library(GENESPACE)
library(rtracklayer)
library(RColorBrewer)
library(ggbeeswarm)
library(Gviz)
library(GenomicRanges)
library(viridis)
library(VariantAnnotation)
library(GenomicFeatures)
library(plyranges)
library(circlize)
library(gggenomes)
library(ggnewscale)
library(colorspace)
library(cowplot)
library(patchwork)
source('functions.R')
binwidth = 1e6

################################################################################

## get genome sizes

################################################################################

fai <- fread('genespace/fai/schLug1.fa.fai')
sum(fai$V2)

fais <- list.files('genespace/fai', pattern = '*.fa.fai')

gsize <- sapply(paste0('genespace/fai/',fais), genome_size_from_fai) |>
  as.data.frame() |>
  rownames_to_column(var = 'genome') |>
  mutate(genome = gsub('genespace/fai/(sch.+).fa.fai','\\1',genome))

names(gsize) <- c('genome', 'size')


################################################################################
## The matching blocks between genomes
genome_names <- c(
  'schMedS3h1',
  'schMedS3h2',
  'schPol2',
  'schLug1',
  'schNov1'
)

block_coords <- fread('genespace/results/blkCoords.txt')
phased_blocks <- fread('genespace_/riparian/refPhasedBlkCoords.txt')

gffs <- list(
  'EDTA.anno.bed/schMedS3_h1.fa.mod.EDTA.TEanno.gff3.gz',
  'EDTA.anno.bed/schMedS3_h2.fa.mod.EDTA.TEanno.gff3.gz',
  'EDTA.anno.bed/schPol2.fa.mod.EDTA.TEanno.gff3.gz',
  'EDTA.anno.bed/schLug1.fa.mod.EDTA.TEanno.gff3.gz',
  'EDTA.anno.bed/schNov1.fa.mod.EDTA.TEanno.gff3.gz'
)

GFFS <- lapply(gffs, import)
DFS <- lapply(GFFS, as.data.frame)

for (i in seq_along(DFS)) {
  DFS[[i]]$genome <- genome_names[i]
}

DFS[[5]] <- DFS[[5]] |>
  mutate(seqnames = gsub('^([0-9])$', 'chr_\\1', seqnames))
joint_df <- do.call(rbind, DFS)

rm(DFS)

# remove the satDNA annotation since we have a newer custom version

joint_df <- joint_df[!grepl('CL', joint_df$Name),]

################################################################################

## Also load our satellite annotation

################################################################################

gffs <- list(
  '02_combined_satDNA/schMedS3_h1_satDNA_min1k_L-200.gff3',
  '02_combined_satDNA/schMedS3_h2_satDNA_min1k_L-200.gff3',
  '02_combined_satDNA/schPol2_satDNA_min1k_L.gff3',
  '02_combined_satDNA/schLug1_satDNA_min1k_L.gff3',
  '02_combined_satDNA/schNov1_satDNA_min1k_L.gff3'
)

GFFS <- lapply(gffs, import)
DFS <- lapply(GFFS, as.data.frame)

for (i in seq_along(DFS)) {
  DFS[[i]]$genome <- genome_names[i]
}

DFS[[5]] <- DFS[[5]] |>
  mutate(seqnames = gsub('^([0-9])$', 'chr_\\1', seqnames))


sat_df <- do.call(rbind, DFS)
sat_df$Classification = 'satDNA'
rm(DFS)

################################################################################
## join the in a simple dataframe for plotting
sat_plt_df <- select(sat_df, - original_names, - strands)

plot_df <- joint_df |> 
  select(names(sat_plt_df)) |> 
  rbind(sat_plt_df)

################################################################################

## Simplify names & setup color

################################################################################

plot_df <- plot_df |>
  mutate(simple_group = case_when(
    grepl("hAT|NHAT", Name) ~ "hAT-family",
    grepl("Mariner|MARINER|MAR", Name) ~ "Mariner-family",
    grepl("MERLIN|Merlin", Name) ~ "MERLIN-family",
    grepl("rnd", Name) ~ "rnd-family",
    grepl("SMAR", Name) ~ "SMAR-family",
    grepl("[Pp]iggyBac|PIGGYB", Name) ~ "piggyBac-family",
    grepl("MuDR", Name) ~ "MuDR-family",
    grepl("SLF", Name) ~ "SLF-family",
    grepl("^DNA", Name) ~ "DNA-family",
    grepl("BURRO", Name) ~ "BURRO-family",
    grepl("Polinton|POLINTN", Name) ~ "Polinton-family",
    grepl("TE_00", Name) ~ "EDTA_novel",
    grepl("Satel", Name) ~ "satDNA",
    grepl("chr", Name) ~ "EDTA_chr",
    grepl("TTAGGG", Name) ~ 'Telomere',
    TRUE ~ Name),
    simple_class = case_when(
      grepl("^DNA", Classification) ~ "DNA-element",
      grepl("^TIR", Classification) ~ "DNA-element",
      grepl("^MITE", Classification) ~ "DNA-element",
      grepl("^polin", Classification) ~ "DNA-element",
      grepl("^LINE", Classification) ~ "LINE-element",
      grepl("^LTR", Classification) ~ "LTR-element",
      grepl("^SINE", Classification) ~ "SINE-element",
      simple_group == "SINEL_SM" ~ "SINE-element",
      simple_group == "Polinton-family" ~ "DNA-element",
      simple_group == "PIVE" ~ "DNA-element",
      TRUE ~ simple_group
    )
  )

################################################################################

## TE color scale

################################################################################

element_order <- c(
  "DNA-element" , 
  "LINE-element" ,
  "SINE-element",
  "LTR-element",
  "EDTA_novel",
  "satDNA",
  "Telomere"
)

# colors
element_color <- c(
  "#E7298A",
  "#1B9E77","#D95F02",
  "#9E9AC8", "#7570B3",
  "#66A61E", 'black'
)

names(element_color) <-  element_order

################################################################################

## Genome order

################################################################################

gorder <- c(
  'schMedS3h1',
  'schMedS3h2',
  'schPol2',
  'schNov1',
  'schLug1'
)

################################################################################

## Summarise by 

################################################################################
joint_gsum <- plot_df |>
  group_by(genome,type) |>
  summarise(N = n(),
            total_length = sum(width))

joint_gclass <- plot_df |>
  group_by(genome, Classification) |>
  summarise(N = n(),
            total_length = sum(width))

joint_gsclass <- plot_df |>
  group_by(genome, simple_class) |>
  summarise(N = n(),
            total_length = sum(width))

################################################################################
## add genome size information
joint_gsum$gsize <- gsize$size[match(joint_gsum$genome, gsize$genome)]
joint_gsum$genome_coverage <- joint_gsum$total_length / joint_gsum$gsize 

joint_gclass$gsize <- gsize$size[match(joint_gclass$genome, gsize$genome)]
joint_gclass$genome_coverage <- joint_gclass$total_length / joint_gclass$gsize 

joint_gsclass$gsize <- gsize$size[match(joint_gsclass$genome, gsize$genome)]
joint_gsclass$genome_coverage <- joint_gsclass$total_length / joint_gsclass$gsize

################################################################################
joint_gsclass |>
  ggplot(aes(x = genome, y = genome_coverage, fill=factor(simple_class, levels = names(element_color)))) +
  geom_bar(position="stack",stat = "identity") + 
  scale_fill_manual(values = element_color) + 
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')


# absolute length
unan <- joint_gsclass |>
  group_by(genome, gsize) |>
  summarise(annotated = sum(total_length)) |>
  mutate(total_length = gsize - annotated,
         genome_coverage = total_length/gsize,
         N = 1, simple_class = 'unannotated') |>
  select(genome,  simple_class,       N, total_length,      gsize, genome_coverage)


element_color <- c(element_color, 'grey')
names(element_color)[length(element_color)] <- 'unannotated'


pdf('fig/EDTA_TE_content.pdf')
rbind(joint_gsclass, unan) |>
  filter(genome %in% gorder) |> 
  #filter(simple_class != 'satDNA') |> 
  ggplot(aes(x = factor(genome, levels = gorder), y = total_length/1e6, fill=factor(simple_class, levels = names(element_color)))) +
  geom_bar(position="stack",stat = "identity") + 
  scale_fill_manual(values = element_color) + 
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.title.x = element_blank()) +
  ylab('size in Mb') + 
  ggtitle('Repeat annotation using EDTA & RepeatExplorer')
dev.off()


sum_TE_coverage <- joint_gsclass |> 
  group_by(genome) |>
  summarise(rep_coverage = sum(genome_coverage))


################################################################################
## some sorting
# replace one element
new_order <- c("DNA", "DNA/CMC-Chapaev-3", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/hAT", "DNA/hAT-Blackjack", 
            "DNA/hAT-Charlie", "DNA/hAT-hAT19", "DNA/hAT-hATm", "DNA/hAT-hATx", "DNA/hAT-Tag1", "DNA/hAT-Tip100", "DNA/Helitron",
            "DNA/Maverick", "DNA/Merlin", "DNA/MULE-MuDR", "DNA/PiggyBac", "DNA/TcMar-Fot1", "DNA/TcMar-ISRm11", "DNA/TcMar-Mariner",
            "DNA/TcMar-Pogo", "DNA/TcMar-Stowaway", "DNA/TcMar-Tc1", "DNA/TcMar-Tc2", "DNA/TcMar-Tigger", "DNA/Zator",
            "LINE/Penelope", "LINE/R2", "LINE/unknown",
            "LTR/Copia", "LTR/Gypsy", "LTR/unknown",
            "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", "MITE/DTT",
            "TIR/Merlin", "TIR/PiggyBac", "TIR/Tc1_Mariner",
            "DIRS", 
            "Other/DNA_virus", "Penelope", "polinton", "satDNA",
            "Unknown")

okabe_ito_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

degree = 8
# Create lighter/darker versions of the Okabe-Ito colors
lighter_colors <- sapply(okabe_ito_colors, function(color) colorRampPalette(c(color, "#FFFFFF"))(degree)[degree-1])
darker_colors <- sapply(okabe_ito_colors, function(color) colorRampPalette(c(color, "#000000"))(degree)[degree-4])

image(matrix(1:length(okabe_ito_colors), nrow = 1), col = okabe_ito_colors , xaxt = "n", yaxt = "n")
image(matrix(1:length(lighter_colors ), nrow = 1), col = lighter_colors , xaxt = "n", yaxt = "n")
image(matrix(1:length(darker_colors ), nrow = 1), col = darker_colors , xaxt = "n", yaxt = "n")

# Generate the color vectors
color_vector_DNA <- colorRampPalette(c(okabe_ito_colors[8], lighter_colors[8]))(29) # Orange to lighter Orange
color_vector_LINE <- colorRampPalette(c(okabe_ito_colors[4], lighter_colors[4]))(3) # Bluish Green to lighter Bluish Green
color_vector_LTR <- colorRampPalette(c(okabe_ito_colors[6], lighter_colors[6]))(3) # Blue to lighter Blue
color_vector_MITE <- colorRampPalette(c(okabe_ito_colors[2], darker_colors[2]))(5) # Reddish Purple to darker Reddish Purple
color_vector_TIR <- colorRampPalette(c(okabe_ito_colors[7], darker_colors[7]))(3) # Black to darker Black

length(new_order) - length(color_vector_DNA) - length(color_vector_LINE) - length(color_vector_LTR) - length(color_vector_MITE) - length(color_vector_TIR)

color_vector <- c(color_vector_DNA, color_vector_LINE, color_vector_LTR, color_vector_MITE, color_vector_TIR, 
                  "steelblue",
                  okabe_ito_colors[5], "#FCFBE3", # lighter_colors[5], 
                  "#D41159", "grey30")

################################################################################

## for the unsimplified classes
sort(unique(joint_gclass$Classification))

plot_df2 <- rbind(joint_gclass, unan) |>
  filter(genome %in% gorder) |> 
  mutate(Classification = ifelse(Classification == "Unspecified", "Unknown", Classification)) 


pdf('fig/EDTA_TE_content_Classification.pdf', width = 10, height = 10)
plot_df2 |> 
  ggplot(aes(x = factor(genome, levels = gorder), y = total_length/1e6, fill=factor(Classification, levels = new_order))) +
  geom_bar(position="stack",stat = "identity", width = 0.5) + 
  scale_fill_manual(values = color_vector) + 
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.title.x = element_blank()) +
  ylab('size in Mb') + 
  guides(color = guide_legend(ncol = 1)) +
  ggtitle('Repeat annotation using EDTA & RepeatExplorer') +
  theme(legend.position = 'right',
        legend.text = element_text(size = 12),  # Change text size to 12
        legend.title = element_text(size = 14))  # Change title size to 14
dev.off()

pdf('fig/EDTA_TE_content_Classification_nolegend.pdf', width = 3, height = 5)
plot_df2 |> 
  ggplot(aes(x = factor(genome, levels = gorder), y = total_length/1e6, fill=factor(Classification, levels = new_order))) +
  geom_bar(position="stack",stat = "identity", width = 0.5) + 
  scale_fill_manual(values = color_vector) + 
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.title.x = element_blank()) +
  ylab('Genome size (Mb)') + 
  guides(color = guide_legend(ncol = 1)) +
  ggtitle('Repeat annotation using EDTA & RepeatExplorer') +
  theme(legend.position = 'none',
        legend.text = element_text(size = 12),  # Change text size to 12
        legend.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  # Change title size to 14
dev.off()



################################################################################

## Strata

################################################################################
S3_strata <- read.table('../2021-hic/strata/S3h1_strata_from_A2h1_cactus.tsv')
cl <- RColorBrewer::brewer.pal(n = 11, name = "Set3")
S3_strata$colors <- c(cl, cl[11])
################################################################################

## setup for density calculations

################################################################################
seqs <- gggenomes::read_fai("genespace/fai//schMedS3h1.fa.fai") %>%
  filter(grepl("chr", seq_id))
glen <- seqs$length
names(glen) <- seqs$seq_id
################################################################################

## setup parameters

################################################################################
chrcol <- RColorBrewer::brewer.pal(4, "Dark2")
names(chrcol) <- paste0('chr',1:4)
accent_col <- c("#a0cbd8", "#785F59", "#fadb6a", "#a7ada6")
accent_alpha=1
border_col = 'white'
lighten(chrcol)
pal = c(chrcol, 
        rev(lighten(chrcol, amount = 0.5)))

col_text <- "black"
col_chr <- "steelblue"
synal_col = "grey30"
sml_height = 0.05
big_height = 0.3
axis_u <- 25e6
axis_breaks <- seq(0, max(seqs$length) + axis_u, axis_u)

chromosomes <- seqs |>
  dplyr::rename(chr = seq_id, end = length) |>
  mutate(start = 0) |>
  select(chr, start, end)


################################################################################

## overlap of gaps with repeats

################################################################################
TE_S3h1 <- plot_df |> filter(genome == 'schMedS3h1')
# TE_S3h1_chr <- TE_S3h1[grepl('chr',seqnames(TE_S3h1)),]
################################################################################

## plotting

################################################################################
# Define a custom function to process and plot genomic tracks

hex_values <- c("#feda75", "#fa7e1e", "#d62976", "#4f5bd5")
hex_values <- c("#E0E0E0", "#BDBDBD", "#969696", "#737373")
hex_values <- accent_col


chrcol2 <- hex_values
names(chrcol2) <- c('chr1', 'chr2', 'chr3','chr4')

process_genomic_track <- function(block_coords, genome1, genome2) {
  print(genome1)
  print(genome2)
  data <- block_coords %>%
    filter(genome1 == !!genome1,
           genome2 == !!genome2) %>%
    select(chr1, startBp1, endBp1, everything()) %>%
    unique()
  
  gr <- makeGRangesFromDataFrame(data, seqnames.field = 'chr1', start.field = 'startBp1', end.field = 'endBp1',
                                  keep.extra.columns = TRUE)
  
  cl <- chrcol2[gsub('_h2', '', data$chr2)]
  data$color <- ifelse(is.na(cl), adjust_transparency(accent_col[[1]], alpha = 1), cl)
  #data$color <- "#a7ada6"
  circos.genomicTrack(
    data, stack = TRUE, track.height = 0.07, bg.border=NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = value$color, border = border_col)
    })
  
  
  return(gr)
}

# Define genome pairs
genome_pairs <- list(c("schMedS3h1", "schMedS3h2"),
                     c("schMedS3h1", "schPol2"),
                     c("schMedS3h1", "schNov1"),
                     c("schMedS3h1", "schLug1"))


#pdf('fig/genespace_gaps_circo_S3h1.pdf', width = 14, height = 14)
pdf(file = 'fig/preakpoints_TEs.pdf', width = 14, height = 14)
circos.clear()
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.02), start.degree = 85, gap.degree =c(2,2,2,10))
circos.genomicInitialize(chromosomes, plotType = NULL)

# genomes
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  
  # axis text can remain
  circos.text(mean(xlim), mean(ylim), chr, cex=1, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
  
  circos.axis(h="top", 
              labels.cex=0.8, 
              major.at=axis_breaks, 
              labels=round(axis_breaks/1e6, 1),
              col="black", 
              labels.col="black", 
              lwd=0.8, 
              labels.facing="clockwise")
}, bg.col=pal, bg.border=NA, track.height=0.06)

# Process and plot genomic tracks for each genome pair
pair = genome_pairs[[1]]
gr_list <- lapply(genome_pairs, function(pair) {
  process_genomic_track(block_coords, pair[1], pair[2])
})
dev.off()

# Assign individual GRanges objects to separate variables
gr_S3h2 <- gr_list[[1]]
gr_Spol <- gr_list[[2]]
gr_Snov <- gr_list[[3]]
gr_Slug <- gr_list[[4]]

################################################################################

## What is in the gaps?

################################################################################

## Zoom on the gaps

################################################################################

################################################################################
## check what is in the breaks on schmidtea
# get TE annotation
TE_S3h1_chr <- plot_df |> 
  filter(genome == 'schMedS3h1', grepl('chr', seqnames)) |> 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

################################################################################
## get the gaps in S3h2
gap_df <- gaps(gr_S3h2) |>
  as.data.frame()

# make a subset df that has a name to link it with plotting data and internally the same coordinates
gap_df$gapID <- paste(gap_df$seqnames,gap_df$start,gap_df$end, sep = '-')

################################################################################
## TEs in the gaps
gap_TE <- subsetByOverlaps(TE_S3h1_chr, gaps(gr_S3h2))

# overlap for each gap
ol <- findOverlaps(query = TE_S3h1_chr, subject = gaps(gr_S3h2))

gap_TE$color <- element_color[match(gap_TE$simple_class, names(element_color))]
# use overlap info to assign the zoom names
gap_TE$gapID <- gap_df$gapID[subjectHits(ol)]


gap_coverage <- function(TE, gr, mean_coverage) {
  ################################################################################
  
  ## Zoom on the gaps
  
  ################################################################################
  gap_df <- gaps(gr) |>
    as.data.frame()
  # make a subset df that has a name to link it with plotting data and internally the same coordinates
  gap_df$gapID <- paste(gap_df$seqnames,gap_df$start,gap_df$end, sep = '-')
  ################################################################################
  ## TEs in the gaps
  gap_TE <- subsetByOverlaps(TE, gaps(gr))
  # overlap for each gap
  ol <- findOverlaps(query = TE, subject = gaps(gr))
  gap_TE$color <- element_color[match(gap_TE$simple_class, names(element_color))]
  # use overlap info to assign the zoom names
  gap_TE$gapID <- gap_df$gapID[subjectHits(ol)]
  
  
  sum_class <- gap_TE |>
    as.data.frame() |>
    group_by(gapID, simple_class) |>
    summarise(N = n(),
              total_length = sum(width))
  
  gap_df <- left_join(gap_df, sum_class, by = 'gapID') |>
    mutate(gap_coverage = total_length/width)
  
  gap_df |>
    group_by(gapID) |>
    summarise(pro = sum(gap_coverage))
  
  gap_df |> 
    filter(simple_class == "satDNA", gap_coverage > 0.5)
  
  pp <- gap_df |>
    filter(grepl('[Cc]hr', seqnames)) |>
    ggplot(aes(x = gapID, y = gap_coverage, fill=factor(simple_class, levels = names(element_color)))) +
    geom_bar(position="stack",stat = "identity") + 
    geom_hline(yintercept = mean_coverage) +
    scale_fill_manual(values = element_color) + 
    facet_wrap(~ seqnames, scales = 'free_x') +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = 'bottom',
          axis.text.x = element_blank())
  return(pp)
}

names(GFFs)
################################################################################

## gap coverage in relation to other genomes

################################################################################

co <- sum_TE_coverage[4,] |> pull(rep_coverage)

gap_coverage(TE_S3h1_chr, gr_S3h2, co)
gap_coverage(TE_S3h1_chr, gr_A2h1, co)
gap_coverage(TE_S3h1_chr, gr_A2h2, co)

gap_coverage(TE_S3h1_chr, gr_Spol, co)
gap_coverage(TE_S3h1_chr, gr_Slug, co)
gap_coverage(TE_S3h1_chr, gr_Snov, co)


################################################################################
## there are some gaps that are entierly satDNA 
## we will highlight those.

gr <- makeGRangesFromDataFrame(gap_df, seqnames.field = 'seqnames', start.field = 'start', end.field = 'end',
                               keep.extra.columns = TRUE)

data = gap_df |> 
  filter(simple_class == "satDNA", gap_coverage > 0.5)

cl <- chrcol2[gsub('_h2', '', data$chr2)]
data$color <- 'black'
#data$color <- "#a7ada6"
circos.genomicTrack(
  data, stack = TRUE, track.height = 0.07, bg.border=NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = value$color, border = border_col)
  })



grs <- list(gr_S3h2, gr_Spol, gr_Slug, gr_Snov)
gr_names <- list('S3h2', 'Spol', 'Slug', 'Snov')
current_gr <- gr_Spol
current_name = 'Spol'

refs <- list('schMedS3h1', 'schMedS3h2', 'schPol2', 'schNov1', 'schLug1')
ref = 'S3h1'

TE_S3h1_chr <- plot_df |> 
  filter(genome == 'schMedS3h1', grepl('chr', seqnames)) |> 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)



make_subset_gr(block_coords, 'schNov1', 'schPol2')

################################################################################

## What is immediately flanking in Smed

################################################################################
# Load the necessary libraries
library(future)
library(furrr)

# modify plot_df for schPol2

plot_df3 <- plot_df |> 
  mutate(seqnames = case_when(
    genome == 'schPol2' ~ gsub('^chr_', 'chr', seqnames),
    TRUE ~ seqnames
  ))

# Set up parallelization plan
plan(multisession, workers = 5)
options(future.globals.maxSize= 8912896000)
w = 1e4
nsim = 1000
################################################################################
## we test for each pair for enrichment for TEs in both genomes.
## we are focusing on comparisons that include schMedS3h1

genome_pairs <- list(
  c("schMedS3h1", "schMedS3h2"),
  c("schMedS3h2", "schMedS3h1"),
  
  c("schMedS3h1", "schNov1"),
  c("schNov1", "schMedS3h1"),
  
  c("schMedS3h1", "schLug1"),
  c("schLug1","schMedS3h1"),
  
  c("schMedS3h1", "schPol2"),
  c("schPol2", "schMedS3h1"))

make_subset_gr <- function(block_coords, g1, g2) {
  # to avoid duplication not every genome is in genome1 column
  # this function returns a GRranges object of the blocks in genome1 one
  # check if it is present in the default orientation
  print(paste0('making GRange for comparison of ', g1, ' and ', g2))

  data <- block_coords %>%
    filter(genome1 == !!g1,
           genome2 == !!g2) %>%
    select(chr1, startBp1, endBp1, everything()) %>%
    unique()
  
  data2 <- block_coords %>%
    filter(genome1 == !!g2,
           genome2 == !!g1) %>%
    select(chr2, minBp2, maxBp2, everything()) %>%
    unique()
  
  if (nrow(data) > 0 & nrow(data2) == 0) {
    print('present in default orientation')
    gr <- makeGRangesFromDataFrame(data, seqnames.field = 'chr1', start.field = 'startBp1', end.field = 'endBp1',
                                   keep.extra.columns = TRUE)
  } 
  
  if (nrow(data2) > 0 & nrow(data) == 0) {
    print('using alternative orientation')
    gr <- makeGRangesFromDataFrame(data2, seqnames.field = 'chr2', start.field = 'minBp2', end.field = 'maxBp2',
                                   keep.extra.columns = TRUE)
  } 
  return(gr)
}

simulate_windows <- function(x, w, up, down, chr_to_sample, chromsize, TEs)  {
  
  # sample function
  get_random_starts <- function(seqname, obs, w = 1e3, chromsize) { 
    # number of observations for that seq
    observations = sum(seqnames(obs) == seqname)
    sample(seq(1, chromsize[seqname] - w + 1), observations)
  }
  
  obs = up
  # get starting points
  random_starts_up <- unlist(
    lapply(chr_to_sample, get_random_starts, obs = obs, w = 1e3, chromsize = chromsize)
  )
  
  # make the start points into random intervals
  random_up <- GRanges(
    seqnames = rep(chr_to_sample, times = sapply(chr_to_sample, function(seqname) sum(seqnames(obs) == seqname))),
    ranges = IRanges(start = random_starts_up, width = w)
  )
  
  obs = down
  # get starting points
  random_starts_down <- unlist(
    lapply(chr_to_sample, get_random_starts, obs = obs, w = 1e3, chromsize = chromsize)
  )
  
  # make the start points into random intervals
  random_down <- GRanges(
    seqnames = rep(chr_to_sample, times = sapply(chr_to_sample, function(seqname) sum(seqnames(obs) == seqname))),
    ranges = IRanges(start = random_starts_down, width = w)
  )
  
  # Combine the random intervals.
  random_flanking <- c(random_down, random_up)
  
  # Now you can intersect these random intervals with your TE data to generate an expected count.
  random_flanking_TE <- join_overlap_intersect_within(TEs, random_flanking)
  
  te_random_flanking <- data.frame(random_flanking_TE)
  table_sim <- table(te_random_flanking$Classification)
  
  return(table_sim)
}

for (p in genome_pairs) {
  reference_genome = p[[1]]
  target_genome = p[[2]]
  print(paste0(reference_genome, ' vs ', target_genome))

  ################################################################################
  ## break points in the reference
  current_gr <- make_subset_gr(block_coords, reference_genome, target_genome)
  ref = reference_genome
  current_name = target_genome
  ################################################################################
  ## TE content of reference
  TEs <-  plot_df3 |> 
    filter(genome == reference_genome) |> 
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  ################################################################################
  ## get chromsize
  sizes <- read.table(paste0(reference_genome, '.chromsize'))
  chromsize <- setNames(sizes$V2, sizes$V1)
  ################################################################################
  ## make flanking windows and get their content
  down <- flank(current_gr, width = w, start=TRUE)
  down$flank_width <- w
  down$gapID <- paste('down', seqnames(down), start(down), end(down), sep = '-')
  
  up <- flank(current_gr, width = w, start=FALSE)
  up$flank_width <- w
  up$gapID <- paste('up', seqnames(up), start(up), end(up), sep = '-')
  
  flanking <- c(down, up)
  flanking_TE <- join_overlap_intersect_within(TEs, flanking)
  ################################################################################
  ## testing for significance using resampling
  ## chromosomes to check
  chr_to_sample <- levels(seqnames(current_gr))
  #sims <- lapply(1:nsim, simulate_windows)
  sims <- future_map(1:nsim, ~simulate_windows(.x, w, up, down, chr_to_sample, chromsize, TEs),
                     seed = TRUE)
  #sims <- future_map(1:nsim, simulate_windows)
  sims_lst <- lapply(sims, as.data.frame)
  
  sims_df <- sims_lst |> 
    purrr::reduce(full_join, by = 'Var1')
  names(sims_df) <- c('Var', 1:nsim)
  
  # replace NAs
  sims_df[is.na(sims_df)] <- 0
  
  sim_pdf <- sims_df |> 
    pivot_longer(cols = -Var, names_to = 'replicate', values_to = 'number') 
  
  # Observed Frequency in Flanking Regions
  te_flanking <- data.frame(flanking_TE)
  table_obs <- table(te_flanking$Classification) |> 
    as.data.frame()
  
  names(table_obs) <- c('Var', 'number')
  table_obs$replicate = 'observation'
  
  res_df <- sim_pdf |> 
    group_by(Var) |> 
    left_join(table_obs, by = 'Var') |> 
    transmute(sim_number = number.x,
              obs_number = number.y) |> 
    mutate(obs_number = ifelse(is.na(obs_number),0, obs_number)) |> 
    summarise(simulations = dplyr::n(),
              observed = unique(obs_number),
              mean = mean(sim_number),
              SD = sd(sim_number),
              pval_smaller = sum(obs_number > sim_number) / simulations,
              pval_larger =  sum(obs_number < sim_number) / simulations)
  
  # make a test for general enrichment of TEs
  total_observed <- sum(table_obs$number)
  
  sim_pdf2 <- sim_pdf |> 
    group_by(replicate) |> 
    summarise(number = sum(number))
  
  res_df2 <- sim_pdf2 |> 
    mutate(obs_number = total_observed) |> 
    summarise(simulations = dplyr::n(),
              observed = unique(obs_number),
              mean = mean(number),
              SD = sd(number),
              pval_smaller = sum(obs_number > number) / simulations,
              pval_larger =  sum(obs_number < number) / simulations) |> 
    mutate(Var = 'total')
  
  res_out <- rbind(res_df, res_df2)
  res_out$ref = ref
  res_out$target = current_name
  
  print('making plots')
  pdf(paste0('fig/breakpoint_enrichment_',ref, '_vs_', current_name, '_flankingsize_', w, '.pdf'))
  print(
    sim_pdf |> 
      ggplot(aes(x = number)) +
      geom_histogram() +
      geom_vline(data = table_obs,
                 aes(xintercept = number),
                 color = 'red') +
      facet_wrap( ~ Var, scales = 'free')
  )
  dev.off()
  
  pdf(paste0('fig/breakpoint_enrichment_overall_',ref, '_vs_', current_name, '_flankingsize_', w, '.pdf'))
  print(
    sim_pdf2 |> 
      ggplot(aes(x = number)) +
      geom_histogram() +
      geom_vline(xintercept = total_observed,
                 color = 'red')
  )
  dev.off()
  
  write.table(x = res_out,  file = paste0('fig/breakpoint_enrichment_overall_', ref, '_vs_', current_name, '_flankingsize_', w, '.tsv'),
              sep = '\t', row.names = FALSE)
  
  print(paste0('DONE with ', current_name))
}


################################################################################

## Generate results table

################################################################################

res_files <- list.files('fig/', pattern = '*.tsv', full.names = TRUE)

res_tables <- lapply(res_files, read.table, header = TRUE)

adj_tables <- lapply(res_tables,
       function(df) {
         df2 <- df |> 
           filter(mean > 10)
         
         nrow(df2)
         ps <- p.adjust(p = c(df2$pval_smaller, df2$pval_larger), method = 'fdr')
         
         df2$padj_smaller = ps[1:nrow(df2)]
         df2$padj_larger = ps[(nrow(df2)+1):length(ps)]
         
         df3 <- df2 |> 
           select(Var, padj_smaller, padj_larger) 
         res <- left_join(df, df3) 
         
         return(res)
       })

purrr::reduce(adj_tables, rbind) |> 
  write_xlsx(path= 'fig/breakpoint_enrichment.xlsx')

################################################################################

## Vizualization of enrichment

################################################################################

res_tables_tmp <- readxl::read_xlsx('fig/breakpoint_enrichment.xlsx')

col = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col = c("#E69F00", "#009E73",  "#0072B2", "#CC79A7", 'grey')
sel = c('DNA/DTA', 'DNA/Maverick', 'DNA/PiggyBac', 'LINE/R2', 'LTR/Gypsy', 'satDNA', 'Unknown')
sel = c('total', 'LINE/R2', 'LTR/Gypsy')


res_tables <- res_tables_tmp |> 
  mutate(ref = factor(ref, levels = c('schMedS3h1', 'schMedS3h2', 'schPol2', 'schNov1', 'schLug1')),
  target = factor(target, levels = c('schMedS3h1', 'schMedS3h2', 'schPol2', 'schNov1', 'schLug1'))
)

plot_df <- res_tables |> 
  filter(ref == 'schMedS3h1') |> 
  filter(Var %in% sel) |> 
  mutate(Var = factor(Var, levels = sel)) |> 
  #filter(!is.na(padj_smaller)) |> 
  mutate(deviation = observed - mean,
         sd_low = mean - SD,
         sd_high = mean + SD,
  alpha_value = ifelse(padj_smaller < 0.05 | padj_larger < 0.05, 1, 0.7)) # Adjust alpha based on condition


p1 <- plot_df |> 
  ggplot(aes(x = target, color = target,  alpha = alpha_value)) +
  geom_point(aes(y = mean), stat = 'identity', size = 2) +
  geom_errorbar(aes(ymin = sd_low, ymax = sd_high), width = 0.5) +
  geom_point(aes(y = observed),color = 'black', stat = 'identity', size = 2) +
  scale_color_manual(values = col) +
  scale_alpha_identity() + # Ensures alpha values are used as provided
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_wrap(~Var, nrow = 1, scales = 'free_x') +
  coord_flip()


################################################################################
## For the other species
plot_df2 <- res_tables |> 
  filter(target == 'schMedS3h1') |> 
  filter(Var %in% sel) |> 
  #filter(!is.na(padj_smaller)) |> 
  mutate(Var = factor(Var, levels = sel)) |> 
  mutate(deviation = observed - mean,
         sd_low = mean - SD,
         sd_high = mean + SD,
         alpha_value = ifelse(padj_smaller < 0.05 | padj_larger < 0.05, 1, 0.7)) # Adjust alpha based on condition


p2 <- plot_df2 |>       
  ggplot(aes(x = ref)) +
  geom_point(aes(y = mean, color = ref), size = 2) +
  geom_errorbar(aes(ymin = sd_low, ymax = sd_high, color = ref), width = 0.5) +
  geom_point(aes(y = observed), alpha = 1, size = 2) +
  scale_color_manual(values = col) +
  scale_alpha_identity() + # Ensures alpha values are used as provided
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  facet_wrap(~Var, nrow = 1, scales = 'free_x') +
  coord_flip()


p1 / p2


################################################################################
## combine into one plot with target and ref as facets

plot_df3 <- res_tables |> 
  filter(Var %in% sel) |> 
  mutate(Var = factor(Var, levels = sel)) |> 
  mutate(deviation = observed - mean,
         sd_low = mean - SD,
         sd_high = mean + SD) |> 
  mutate(perspective = case_when(
    ref == 'schMedS3h1' ~ 'schMedS3h1',
    target == 'schMedS3h1' ~ 'target',
    TRUE ~ 'ERROR'
    )) |> 
  mutate(genome = case_when(
    perspective == 'schMedS3h1' ~ target,
    perspective == 'target' ~ ref
    )) |> 
  mutate(
    perspective = factor(perspective, levels = c('schMedS3h1', 'target')),
    genome = factor(genome, levels = rev(c('schMedS3h2', 'schPol2','schNov1', 'schLug1')))
    )


p4 <- plot_df3 |> 
    ggplot(aes(x = genome)) +
    geom_point(aes(y = mean, color = genome), size = 2) +
    geom_errorbar(aes(ymin = sd_low, ymax = sd_high, color = genome), width = 0.5) +
    geom_point(aes(y = observed), alpha = 1, size = 2) +
    scale_color_manual(values = col) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90)) +
    facet_wrap(Var~ perspective, nrow = 1, scales = 'free_x') +
    coord_flip()

ggsave('fig/breakpoint_enrichment_overall_main.svg', p4)
fontS = 6

pT <- plot_df3 |> 
  filter(Var == 'total') |> 
  ggplot(aes(x = genome)) +
  geom_point(aes(y = mean, color = genome), size = 0.2, pch = 3) +
  geom_errorbar(aes(ymin = sd_low, ymax = sd_high, color = genome), width = 0.5) +
  geom_point(aes(y = observed), alpha = 1, size = 0.2, pch = 1) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        text = element_text(size = fontS ),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 2, margin = margin()),
        strip.background = element_rect(
          color="black", size=0.1, linetype="solid"
        )) +
  facet_wrap(Var~ perspective, nrow = 1, scales = 'free_x') +
  coord_flip()
 

pO <- plot_df3 |> 
  filter(Var != 'total') |> 
  ggplot(aes(x = genome)) +
  geom_point(aes(y = mean, color = genome), size = 0.2, pch = 3) +
  geom_errorbar(aes(ymin = sd_low, ymax = sd_high, color = genome), width = 0.5) +
  geom_point(aes(y = observed), alpha = 1, size = 0.2, pch = 1) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        text = element_text(size = fontS ),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 2, margin = margin()),
        strip.background = element_rect(
          color="black", size=0.1, linetype="solid"
        )) +
  facet_wrap(Var~ perspective, nrow = 1, scales = 'free_x') +
  coord_flip()

pB <- pT / pO

ggsave('fig/breakpoint_enrichment_overall_main.pdf', pB,
       width = 70, height = 80, units = 'mm')

library(ggplot2)
library(dplyr)

plot_df <- plot_df |>
  mutate(alpha_value = ifelse(padj_smaller < 0.05 | padj_larger < 0.05, 1, 0.3))

ggplot(plot_df, aes(x = Var, y = mean, fill = target, alpha = alpha_value)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = mean - SD, ymax = mean + SD), position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = okabe_ito_colors) +
  scale_alpha_identity() +
  theme_minimal() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ref)

################################################################################

## There is enrichment for a particular elements

################################################################################

# we should check the interchromosomal breaks

# For this I need to categorize them.
# I will use Spol for testing

# give them a name
gr_Spol$blkID

aste(gap_df$seqnames,gap_df$start,gap_df$end, sep = '-')

grs <- list(gr_S3h2, gr_Spol, gr_Snov, gr_Slug)

sp_names <- unlist(lapply(genome_pairs, '[[', 2))

down_lst = list()

for (i in 1:length(grs)) {
  current_gr <- grs[[i]]
  cat_df <- current_gr |> 
    as.data.frame() |> 
    filter(grepl('chr', seqnames)) |> 
    arrange(seqnames, start) 
    
  for (r in 1:nrow(cat_df)) {
    cat_df[r, 'down_stream_gap'] = case_when(
      cat_df[r, 'seqnames'] != cat_df[r + 1, 'seqnames'] ~ 'query_intra',
      cat_df[r, 'chr2'] == cat_df[r + 1, 'chr2'] ~ 'inter',
      TRUE ~ 'intra'
    )
  }
    
  cat_df <- cat_df |> 
    filter(down_stream_gap != 'query_intra')
  cat_df$species = sp_names[i]
  
  
  down_lst[[i]] <- cat_df
}

all_cat_df <- do.call(rbind, down_lst)


################################################################################
## Now get the flanking content for each category
w = 1e5
gr_all <- c(gr_S3h2, gr_Spol, gr_Snov, gr_Slug)

gr_all$down_stream_gap <- all_cat_df$down_stream_gap[match(gr_all$blkID, all_cat_df$blkID)]

down <- flank(gr_all, width = w, start=TRUE)
down$flank_width <- w
down$gapID <- paste('down', seqnames(down), start(down), end(down), sep = '-')



down_TE <- join_overlap_intersect(down, TE_S3h1_chr)

sp_names

p_gap <- data.frame(down_TE) |> 
  group_by(genome2, blkID, down_stream_gap, simple_class) |>
  filter(simple_class != 'Telomere') |> 
  summarise(N = n(),
            total_length = sum(width),
            fraction = sum(width)/w) |>
  na.omit()  |> 
  group_by(genome2, down_stream_gap, simple_class) |>
  summarise(fraction = median(fraction)) |> 
  ggplot(aes(x = down_stream_gap , y = fraction, fill = simple_class)) +
  geom_bar(position = 'stack',stat = "identity") +
  scale_fill_manual(values = element_color) + 
  scale_x_discrete(labels = c('inter-chr', 'intra-chr')) +
  xlab('') +
  ylab('Coverage') +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
  ) +
  facet_wrap( ~ factor(genome2, levels = sp_names), nrow = 1, labeller = as_labeller(setNames(sp_names, sp_names)))



pdf('fig/breakpoint_TE.pdf')
p_gap
dev.off()



