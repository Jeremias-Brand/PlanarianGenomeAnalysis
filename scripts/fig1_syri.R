library(VariantAnnotation)
library(data.table)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicFeatures)
library(plyranges)
library(circlize)
library(gggenomes)
library(ggnewscale)
library(colorspace)
library(viridis)
binwidth = 1e6


# minimap2 -x asm5 --eqx ref.fa query.fa > minimap2.paf

flip_coordinates <- function(x, xrange) {
  return(xrange[2] - x + xrange[1])
}

split_n_flip <- function(data) {
  region1 <- data[,1:3]
  region2 <- data[,9:11]
  
  chromosome_ranges <- region2 %>%
    group_by(ChrB) %>%
    summarize(min_start = min(StartB), max_end = max(EndB)) 
  
  region2_flipped <- region2 %>%
    left_join(chromosome_ranges, by = "ChrB") %>%
    rowwise() %>%
    mutate(
      StartB = flip_coordinates(StartB, c(min_start, max_end)),
      EndB = flip_coordinates(EndB, c(min_start, max_end))
    ) %>%
    select(-min_start, -max_end) # Remove the extra columns
  
  return(list(region1, region2_flipped))
}


rev_x = function(x, xrange = CELL_META$xlim) {
  xrange[2] - x + xrange[1]
}


calc_den <- function(GR, binwidth = 1e6) {
  cov <- coverage(GR)
  bins <- tileGenome(glen, tilewidth = binwidth, cut.last.tile.in.chrom = TRUE)
  den <- binnedAverage(bins, cov, varname = 'den')
  res <-  as.data.frame(den) %>%
    rename(seq_id = seqnames) %>%
    mutate(strand = ".", scaled = range01(den))
  return(res)
}


makechromsomes <- function(chrom_first, chrom_second, chr_first = "chr", chr_second = chr_first) {
  h1 <- read.table(chrom_first)
  names(h1) <- c("chrom", "size")
  h1 <- filter(h1, grepl(chr_first, chrom)) %>%
    arrange(chrom)
  h1_xlim <- matrix(c(rep(0, nrow(h1)), h1$size), ncol=2)
  
  h2 <- read.table(chrom_second)
  names(h2) <- c("chrom", "size")
  h2 <- filter(h2, grepl(chr_second, chrom)) %>%
    arrange(desc(chrom))
  h2_xlim <- matrix(c(rep(0, nrow(h2)), h2$size), ncol=2)
  
  genome <- rbind(h1,h2)
  genome_xlim <- matrix(c(rep(0, nrow(genome)), genome$size), ncol=2)
  
  chromosomes <- genome %>%
    mutate(chr = chrom, start = 0, end = size) %>%
    select(chr, start, end)
  return(chromosomes)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}



density_ggenomes <- function(seqs, data, binwidth=1e5, varname = "feat_density") {
  # working with gggenome imported files
  # seqs:  gggenomes::read_fai("genome.fa.fai") %>%
  # filter(grepl("Chr", seq_id))
  # data: gggenomes::read_bed(bed) %>% filter(grepl("Chr", seq_id))
  glen <- seqs$length
  names(glen) <- seqs$seq_id
  
  satGR <- data %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seq_id",
                                            seqinfo = glen )
  cov <- coverage(satGR)
  bins <- tileGenome(glen, tilewidth=binwidth, cut.last.tile.in.chrom=TRUE)
  # this object is useful in the GenomicRanges universe but not good for gggenomes
  den <- binnedAverage(bins, cov, varname = varname)
  # We need a few tweaks
  df <- as.data.frame(den) %>%
    dplyr::rename(seq_id = seqnames)
  # only specific strand symbols are allowed.
  df$strand <- "."
  
  return(df)
}


combine_chromosomes <- function(chrom_first, chrom_second,
                       chr_first = "Chr", chr_second = "Chr", rev = TRUE) {
  h1 <- read.table(chrom_first)
  names(h1) <- c("chrom", "size")
  h1 <- filter(h1, grepl(chr_first, chrom)) %>%
    arrange(chrom)
  h1_xlim <- matrix(c(rep(0, nrow(h1)), h1$size), ncol=2)
  
  h2 <- read.table(chrom_second)
  names(h2) <- c("chrom", "size")
  if (rev) {
    h2 <- filter(h2, grepl(chr_second, chrom)) %>%
      arrange(desc(chrom))
  } else {
    h2 <- filter(h2, grepl(chr_second, chrom)) %>%
      arrange(chrom)
  }

  h2_xlim <- matrix(c(rep(0, nrow(h2)), h2$size), ncol=2)
  
  genome <- rbind(h1,h2)
  genome_xlim <- matrix(c(rep(0, nrow(genome)), genome$size), ncol=2)
  
  return(genome)
}


color_by_value <- function(values, n_colors = 25, option = 'magma', trim = 0, trim_col = 'grey', top = FALSE) {
  # trim allows to remove the top or bottom percentiles from the color scale
  # assigning them the trim_col
  # Generate a color palette with n_colors color steps
  my_colors <- viridis(n_colors, option = option, direction = -1)
  # Map the values to colors using the color scale
  my_color_vector <- my_colors[findInterval(values, seq(0, 1, length.out = n_colors + 1), all.inside = TRUE)]
  
  if (trim > 0) {
    # Calculate the trim limits
    lower_limit <- quantile(values, probs = trim)
    upper_limit <- quantile(values, probs = 1 - trim)
    
    # Identify the indices of values within the trimmed regions
    lower_trim_indices <- which(values <= lower_limit)
    upper_trim_indices <- which(values >= upper_limit)
    
    # Replace the corresponding colors with the trim_col
    if (top) {
      my_color_vector[upper_trim_indices] <- trim_col
      within_limits_indices <- which( values < upper_limit)
    } else {
      my_color_vector[lower_trim_indices] <- trim_col
      within_limits_indices <- which(values > lower_limit)
    }
    
    # Map the values to colors using the color scale
    my_color_vector[within_limits_indices] <- my_colors[findInterval(values[within_limits_indices], seq(0, 1, length.out = n_colors + 1), all.inside = TRUE)]
    return(my_color_vector)
  }
  
  return(my_color_vector)
}



simplify_te_ann <- function(df) {
  df2 <- df |>
    mutate(
      simple_group = case_when(
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
        grepl("TE_00", Name) ~ "Unknown",
        grepl("Satel", Name) ~ "satDNA",
        grepl("chr", Name) ~ "Unknown",
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
  return(df2)
}


################################################################################

## coplotting two assemblies in most cases the two phases g1 and g2

################################################################################

g1 = 'schMedS3_h1'
g2 = 'schMedS3_h2'

HCONF = TRUE
################################################################################

## loading the different types of annotation

################################################################################

################################################################################
## For the plotting we need to combine the data on both genomes
## Load genome information & combine

genome = combine_chromosomes(
  chrom_first = "lib/schMedS3_h1.chromsize",
  chrom_second = "lib/schMedS3_h2.chromsize",
  chr_first = "chr", chr_second = "chr")

# rename for plotting
chromosomes <- genome %>%
  mutate(chr = chrom, start = 0, end = size) %>%
  select(chr, start, end)

# named vector of the chromosome lenghts
glen <- chromosomes$end
names(glen) <- chromosomes$chr

################################################################################
## gene annotation

g1_gene <- import.gff3(paste0('latest_annotation/subsets/', g1, '_ENCODE_hybrid_longest_coding_sorted.gff3'))

g1_gene_df <- g1_gene |>
  as.data.frame() |> 
  # important to remove the factor for later subsetting
  mutate(seqnames = as.character(seqnames)) |> 
  mutate(Parent = ifelse(as.character(Parent) == "character(0)", NA, as.character(Parent)) )   


if (HCONF) {
  g1_hconf <- import.gff3(paste0('latest_annotation/subsets/', g1, '_ENCODE_hybrid_agat_only_hconf.gff3'))
  
  g1_hconf_df <- g1_hconf |>
    as.data.frame() |> 
    # important to remove the factor for later subsetting
    mutate(seqnames = as.character(seqnames)) |> 
    mutate(Parent = ifelse(as.character(Parent) == "character(0)", NA, as.character(Parent)))
}


g2_gene <- import.gff3(paste0('latest_annotation/subsets/', g2, '_ENCODE_hybrid_longest_coding_sorted.gff3'))

g2_gene_df <- g2_gene |>
  as.data.frame() |> 
  # important to remove the factor for later subsetting
  mutate(seqnames = as.character(seqnames)) |> 
  mutate(Parent = ifelse(as.character(Parent) == "character(0)", NA, as.character(Parent)) )

if (HCONF) {
  g2_hconf <- import.gff3(paste0('latest_annotation/subsets/', g2, '_ENCODE_hybrid_agat_only_hconf.gff3'))
  
  g2_hconf_df <- g2_hconf |>
    as.data.frame() |> 
    # important to remove the factor for later subsetting
    mutate(seqnames = as.character(seqnames)) |> 
    mutate(Parent = ifelse(as.character(Parent) == "character(0)", NA, as.character(Parent)))
}

## Combine them
g_gene_df <- rbind(g1_gene_df, g2_gene_df)
g_hconf_df <- rbind(g1_hconf_df, g2_hconf_df)

################################################################################
## satDNA
g1_sat <- import.gff3(paste0('02_combined_satDNA/', g1, '_satDNA_min1k_L-200.gff3'))
g1_sat_df <- g1_sat |>
  as.data.frame() |> 
  mutate(Name = gsub('Satellite/', '', Name) )

g2_sat <- import.gff3(paste0('02_combined_satDNA/', g2, '_satDNA_min1k_L-200.gff3'))
g2_sat_df <- g2_sat |>
  as.data.frame() |> 
  mutate(Name = gsub('Satellite/', '', Name) )

g_sat_df <- rbind(g1_sat_df, g2_sat_df)

################################################################################
## TEs
# TODO assess if simplification is needed of the classes used
g1_te <- import.gff3(paste0('EDTA.anno.bed/', g1, '.fa.mod.EDTA.TEanno.gff3.gz'))
# remove the satDNA annotation since we have a newer custom version
g1_te <- g1_te[!grepl('CL', g1_te$Name),]   
# create simplified groups, the simplified groups also contain Telomere annotation
g1_te_df <- g1_te |> as.data.frame() |> 
  simplify_te_ann()

g2_te <- import.gff3(paste0('EDTA.anno.bed/', g2, '.fa.mod.EDTA.TEanno.gff3.gz'))
# remove the satDNA annotation since we have a newer custom version
g2_te <- g2_te[!grepl('CL', g2_te$Name),]   
g2_te_df <- g2_te |> as.data.frame() |> 
  simplify_te_ann()

g_te_df <- rbind(g1_te_df, g2_te_df)
################################################################################

## Load additional data here

################################################################################

################################################################################
## Genespace bed file could be loaded if needed
# some parsing bugs at the moment
# BED <- fread('../2022-orthofinder/genespace_v108/results/combBed.txt')
# gggenomes::read_bed('../2022-orthofinder/genespace_v108/results/combBed.txt', skip = 1)

################################################################################
## Genome alignment
syri_vcf = "14_syri_filter/S3h1_S3h2_chrsyri.vcf.gz"
vcf_file <- syri_vcf
syri <- readVcf(vcf_file)
# there is a hirachy of SVs where they are parents and children that is saved in the 
# parent and ID column
rng <- as.data.frame(rowRanges(syri))
rng$ID <- row.names(rng)

syri_df <- cbind(rng,
                 info(syri)) %>%
  as.data.frame() %>%
  select(seqnames, start, END, ID, REF, ALT, QUAL, FILTER,
         ChrB, StartB, EndB, Parent, VarType, DupType) %>%
  rename(end = END) %>%
  mutate(length = abs(start - end),
         lengthB = abs(StartB - EndB),
         TYPE = unlist(ALT))

SVs <- c("<INVAL>",
         "<DUPVAL>",
         "<TRANSAL>",
         "<SYNAL>",
         "<INVTRAL>",
         "<INVDPAL>"
)

# we could also just plot the Svs without a parent i.e. primary ones ?
synal <- syri_df %>% filter(TYPE %in% SVs,
                            length >= 1e4) %>%
  mutate(colors = case_when(
    TYPE == "<INVAL>" ~  "#7570B3",
    TYPE == "<DUPVAL>" ~  "#1B9E77",
    TYPE == "<TRANSAL>" ~  "#D95F02",
    TYPE == "<SYNAL>" ~ "steelblue",
    TYPE == "<INVTRAL>" ~  "#D95F02",
    TYPE == "<INVDPAL>" ~ "#1B9E77"
  ))


################################################################################
## heterozygosity

names <- c('schMedS3_h1', 'schMedS3_h2')
snp_ranges <- list()
het_ranges <- list()
hom_ranges <- list()

for (i in seq_along(names)) {
  name <- names[i]
  snp_ranges[[i]] <- readRDS(file = paste0("GATK_", name, "/", name, "_snp_range.rds"))
  het_ranges[[i]] <- readRDS(file = paste0("GATK_", name, "/", name, "_het_range.rds"))
  hom_ranges[[i]] <- readRDS(file = paste0("GATK_", name, "/", name, "_hom_range.rds"))
}

snp_range <- do.call(c, snp_ranges)
het_range <- do.call(c, het_ranges)
hom_range <- do.call(c, hom_ranges)




# calculate density
snpden <- coverage(snp_range)
bins <- tileGenome(seqinfo(snp_range), tilewidth=1e6, cut.last.tile.in.chrom=TRUE)
gr <- GRanges(seqnames = seqinfo(snp_range))
windows <- unlist(slidingWindows(gr, width = 1e5, step = 1e4))

snpden_binned <- binnedAverage(bins, snpden, "snpdensity")

het_cov <- coverage(het_range)
het_binned <- binnedAverage(bins = bins, numvar = het_cov, varname = "het_snpden", na.rm = FALSE)
het_windows <- binnedAverage(bins = windows, numvar = het_cov, varname = "het_snpden", na.rm = FALSE)


hom_cov <- coverage(hom_range)
hom_binned <- binnedAverage(bins = bins, numvar = hom_cov, varname = "hom_snpden", na.rm = FALSE)
hom_windows <- binnedAverage(bins = windows, numvar = hom_cov, varname = "hom_snpden", na.rm = FALSE)

all_binned <- left_join(
  data.frame(het_binned),
  data.frame(hom_binned)
)


long <- all_binned %>%
  pivot_longer(cols = c(het_snpden, hom_snpden))

long %>%
  filter(grepl(seqnames, pattern = "chr")) %>%
  ggplot(aes(x = end)) +
  geom_line(aes(y = value, color = name)) +
  facet_wrap( ~ seqnames, scales = "free_x") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")


################################################################################

## Heterozygosity

################################################################################

# this is only run for S3h1... we also need to calcualte it for S3h2


################################################################################
## indel density

################################################################################

## loading genome information

################################################################################


################################################################################

## density calculations

################################################################################

################################################################################
## subset to only chromosomes
chrom = 'chr'
# join data from both chromosomes
g_gene_df <- rbind(g1_gene_df, g2_gene_df)

GR <- g_gene_df %>%
  filter(grepl(chrom, seqnames )) |>
  filter(type == 'gene') |> 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen)

gene_den <- calc_den(GR)
################################################################################
## hconf
GR <- g_hconf_df %>%
  filter(grepl(chrom, seqnames )) |>
  filter(type == 'gene') |> 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen)

hconf_den <- calc_den(GR)

################################################################################
## satDNA
cutoff = 0.3
GR <- g_sat_df %>%
  filter(grepl(chrom, seqnames )) |>
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen)

sat_den <- calc_den(GR)
# binary mask
sat_den$mask <- ifelse(sat_den$scaled >= cutoff, 1, 0)

large_sat <- g_sat_df |> 
  filter(width > 100000) 


################################################################################
## broad TE

GR <- g_te_df %>%
  filter(grepl(chrom, seqnames )) |>
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen)

te_den <- calc_den(GR)

################################################################################
## preparation


################################################################################

## circlize plotting

################################################################################

## Very important to make sure all h2 tracks are reversed! - or nothing is reversed

g1_grep = 'h1'

################################################################################

## pdf

################################################################################

# "#cc79a7" H3K4me3
# 
# "#0072b2" H3K27ac
# 
# "#e69f00" ATAC

col1 <- "#D55E00"
col2 <- "#F0E442"
col3 <- "#009E73"
col4 <- "#56B4E9"
col5 <- "grey70"

################################################################################


helper <- function() {
  circos.clear()
  chrcol <- RColorBrewer::brewer.pal(4, "Dark2")
  lighten(chrcol)
  pal = c(chrcol, 
          rev(lighten(chrcol, amount = 0.5)))
  
  col_text <- "black"
  col_chr <- "steelblue"
  synal_col = "grey30"
  sml_height = 0.05
  big_height = 0.3
  lw = 0.1
  circos.par("track.height"=0.1, 
             #gap.degree=5, 
             #cell.padding=c(0.0, 0, 0, 0.0),
             cell.padding=c(0.002, 0.01, 0.002, 0.001),
             start.degree = 80 ,
             clock.wise = TRUE,
             gap.degree = c(2,2,2,4, 2,2,2,20),
             track.margin = c(0.01,0.01)
             
  )
  
  axis_u <- 20e6
  axis_breaks <- seq(0, ceiling(max(genome$size)/1e6)*1e6, axis_u)
  
  # Function to generate sequence
  generate_sequence <- function(size) {
    breaks <- seq(0, size, by=axis_u )
    print(breaks)
    if (size - breaks[length(breaks)] > 10e6) {
      breaks <- c(breaks, size)
    } else {
      breaks[length(breaks)] <- size
    }
    breaks
  }
  
  # Apply the function to each size
  chromosome_params <- lapply(genome$size, generate_sequence)
  
  # Convert to a list with chromosome names
  names(chromosome_params) <- genome$chrom
  
  ################################################################################
  ## plot chromosomes
  
  circos.genomicInitialize(chromosomes, plotType = NULL)
  # genomes
  circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    
    # axis text can remain
    # circos.text(mean(xlim), mean(ylim), chr, cex=1, col=col_text, 
    #             facing="bending.inside", niceFacing=TRUE)
    params = chromosome_params[[chr]]
    
    seq_include_last <- function(from, to, by) {
      s <- seq(from, to, by)
      if (tail(s, n = 1) != to) {
        s <- c(s, to)
      }
      s
    }
    
    # don't reveres g1
    if(grepl(g1_grep, CELL_META$sector.index)) {
      circos.axis(h="top", 
                  labels.cex=0.8, 
                  major.at=params, 
                  labels=round(params/1e6, 0),
                  col="black", 
                  labels.col="black", 
                  lwd=0.8, 
                  labels.facing="clockwise")
      # reverse g2
    } else {
      circos.axis(h="top", 
                  labels.cex=0.8, 
                  major.at=rev_x(params), 
                  labels=round(params/1e6, 0),
                  col="black", 
                  labels.col="black", 
                  lwd=0.8, 
                  labels.facing="clockwise")
    }
  }, bg.border=F, track.height=0.01)
  
  ################################################################################
  ## Het density
  data <- long %>%
    filter(name == "het_snpden",
           grepl(seqnames, pattern = "chr")) %>%
    mutate(chr = as.character(seqnames))
  
  
  data <- data[,c(1:3, 7)]    
  data$value <- range01(data$value)
  circos.genomicTrack(
    data, track.height = 0.1, track.index = 2,
    ylim=c(0,1),
    panel.fun = function(region, value, ...) {
      if(grepl(g1_grep, CELL_META$sector.index)) {
        circos.genomicLines(region, value, 
                            col = col4,
                            border = col4,
                            area = TRUE,
                            lwd = lw
        )
      } else {
        circos.genomicLines(rev_x(region), value,
                            col = col4,
                            border = col4,
                            area = TRUE,
                            lwd = lw)
      }
    },bg.col = "grey95" )
  
  ################################################################################
  ## HCONF
  data <- hconf_den[,c(1:3,7)]    
  data$colors <- color_by_value(data$scaled, n_colors = 25, option = 'magma', trim = 0.25, trim_col = 'white')
  # Create color scale
  color_scale <- colorRampPalette(c("white", "black"))
  
  # Assign colors
  data$colors <- sapply(data$scaled, function(x) color_scale(100)[ceiling(x * 99) + 1])
  
  
  circos.genomicTrack(
    data, track.height = 0.1,track.index = 3,
    ylim=c(0,1),
    panel.fun = function(region, value, ...) {
      if(grepl(g1_grep, CELL_META$sector.index)) {
        circos.genomicLines(region, value, col = col1, border = col1, area = TRUE, lwd = lw)
      } else {
        circos.genomicLines(rev_x(region), value,  col = col1, border = col1, area = TRUE, lwd = lw)
      }
    },bg.col = "grey95")
  
  ## TE
  data <- te_den[,c(1:3,7)]    
  data$colors <- color_by_value(data$scaled, n_colors = 25, option = 'magma', trim = 0.25, trim_col = 'white')
  
  data$colors <- sapply(data$scaled, function(x) color_scale(100)[ceiling(x * 99) + 1])
  
  
  circos.genomicTrack(
    data, track.height = 0.1,track.index = 4,
    ylim=c(0,1),
    panel.fun = function(region, value, ...) {
      if(grepl(g1_grep, CELL_META$sector.index)) {
        circos.genomicLines(region, value, col = col3 , border = col3 , area = TRUE, lwd = lw)
      } else {
        circos.genomicLines(rev_x(region), value,  col = col3 , border = col3 , area = TRUE, lwd = lw)
      }
    },bg.col = "grey95")
  
}

# bg.col=pal, bg.border=F, track.height=0.02)
###############################

# parents

parents <- syri_df |> 
  filter(Parent == '.') |> 
  mutate(
    colors = case_when(
      TYPE == "<SYN>"    ~ "grey70", 
      TYPE == "<DUP>"    ~ "grey70", 
      TYPE == "<NOTAL>"  ~ color_palette[3],
      TYPE == "<INV>"    ~ col2,
      TYPE == "<INVDP>"  ~ color_palette[5],
      TYPE == "<INVTR>"  ~ col2,
      TYPE == "<TRANS>"  ~ color_palette[7],
      TRUE               ~ color_palette[8] 
    ))  |>
  filter(TYPE != "<NOTAL>",
         TYPE != "<DUP>",
         TYPE != "<INVDP>") 


parents_flip <- split_n_flip(parents)

parents_cols <- scales::alpha(
  parents$colors, alpha=1)



###################

# show the inversion in Chr1

# big inversion
INV81_start <- parents[parents$ID == 'INV81','start']
INV81_end <- parents[parents$ID == 'INV81','end']

align <- syri_df |> 
  # remove small variants
  filter(grepl('<', TYPE)) |> 
  # remove variants that are short
  filter(length > 0.1e6) |> 
  # remove the parents
  filter(Parent != '.') |> 
  mutate(
    colors = case_when(
      TYPE == "<SYNAL>"    ~ "#56B4E9", 
      TYPE == "<DUPAL>"    ~ "grey70", 
      TYPE == "<NOTAL>"  ~ 'red',
      TYPE == "<INVAL>"    ~ "#F0E442",
      TYPE == "<INVDPAL>"  ~ "#F0E442",
      TYPE == "<INVTRAL>"  ~ "#D55E00",
      TYPE == "<HDR>"  ~ "grey30",
      TRUE               ~  'black'
    )) |> 
  filter(seqnames %in% c('chr1_h1'),
         start >= INV81_start, end <= INV81_end,
         (length + lengthB)/ 2 > 5000000) |> 
  filter(TYPE != '<SYNAL>') 
  

align_flip <- split_n_flip(align)

align_cols <- scales::alpha(
  align$colors, alpha=0.8)



pdf('fig/fig_1_v02_het_hconf_te_satDNA_rRNA_syri.pdf', height = 10, width = 10)
helper()

circos.genomicLink(region1 = parents_flip[[1]], region2 = parents_flip[[2]], 
                   lwd = 1,
                   col = parents_cols,
                   border=NA)

circos.genomicLink(region1 = align_flip[[1]], region2 = align_flip[[2]], 
                   lwd = 1,
                   col = "black",
                   border="black")
dev.off()
###################






















################################################################################

## content of insertions in S3h1

################################################################################

df <- read.csv('dotplots/schMedS3_h1__schMedS2_guo_100kb_gaps_in_schMedS3_h1.csv')

large_df <- df |> 
  filter(grepl('chr', seqnames), 
         width >= 1e6)

# overlap with satDNA
# TE
# gene
large_df
g_sat_df


te_GR <- g_te_df %>%
  filter(grepl(chrom, seqnames )) |>
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen,
                                          keep.extra.columns = TRUE)

sat_GR <- g_sat_df %>%
  filter(grepl(chrom, seqnames )) |>
  mutate(Name = gsub('.+/([^/]+$)', '\\1', Original_names)) |> 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen,
                                          keep.extra.columns = TRUE)

srf_GR <- g_srf_df |> 
  filter(grepl(chrom, seqnames )) |>
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen,
                                          keep.extra.columns = TRUE)

G_GR <- large_df %>%
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames" ,
                                          seqinfo = glen)

query = sat_GR
subject = G_GR

find_overlaps <- function(query, subject, type = 'within') {
  overlaps <- findOverlaps(query = query, subject = subject, type = type)
  
  # Get indices of overlapping ranges in query and subject
  query_idx <- queryHits(overlaps)
  subject_idx <- subjectHits(overlaps)
  
  # Extract overlapping ranges from query and subject
  query_overlaps <- query[query_idx]
  subject_overlaps <- subject[subject_idx]
  
  # make subject name
  subject_df <- subject_overlaps |> 
    as.data.frame() |> 
    unite(subject_name, seqnames, start, end, sep = '_', remove = FALSE) |> 
    dplyr::rename(subject_width = width,
           subject_seqnames = seqnames,
           subject_start = start,
           subject_end = end ) |> 
    select(-strand)
  
  query_df <- query_overlaps |> 
    as.data.frame() 
    
    
  overlap_df <- cbind(subject_df, query_df)
  
  return(overlap_df)
}

overlap_df <- find_overlaps(sat_GR, G_GR)

sum_sat <- overlap_df |> 
  group_by(subject_name, subject_width, Name) |> 
  summarise(N=n(), sum_width = sum(width)) |> 
  mutate(coverage = sum_width / subject_width) 

overlap_srf <- find_overlaps(srf_GR, G_GR)

sum_srf <- overlap_srf  |> 
  group_by(subject_name, subject_width, name) |> 
  summarise(N=n(), sum_width = sum(width)) |> 
  mutate(coverage = sum_width / subject_width) 

sum_srf |> 
  ggplot(aes(x = subject_name, y = coverage, fill = name)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE)

ggsave('fig/largest_gaps_srf.pdf')

overlap_te <- find_overlaps(te_GR, G_GR)

sum_te <- overlap_te  |> 
  filter(type != 'Unknown') |> 
  group_by(subject_name, subject_width, type) |> 
  summarise(N=n(), sum_width = sum(width)) |> 
  mutate(coverage = sum_width / subject_width) |> 
  dplyr::rename(Name = type)

s <- rbind(sum_sat, sum_te)

s |> 
ggplot(aes(x = subject_name, y = coverage, fill = Name)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE)

ggsave('fig/largest_gaps_te_satDNA.pdf')


overlap_df |> 
  group_by(subject_name, subject_width) |> 
  summarise(N=n(), sum_width = sum(width)) |> 
  mutate(coverage = sum_width / subject_width) |> 
  ggplot(aes(x = subject_name, y = coverage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal()


overlap_df |> 
  group_by(subject_name, subject_width, Name) |> 
  summarise(N=n(), sum_width = sum(width)) |> 
  mutate(coverage = sum_width / subject_width) |> 
  filter(N > 100) |> 
  ggplot(aes(x = subject_name, y = coverage, fill = Name)) +
  geom_bar(stat = "identity") +
  theme_minimal()

overlap_te <- find_overlaps(te_GR, G_GR)

overlap_te  |> 
  filter(type != 'Unknown') |> 
  group_by(subject_name, subject_width, type) |> 
  summarise(N=n(), sum_width = sum(width)) |> 
  mutate(coverage = sum_width / subject_width) |> 
  filter(N > 100) |> 
  ggplot(aes(x = subject_name, y = coverage, fill = type)) +
  geom_bar(stat = "identity") +
  theme_minimal()

g_sat_df |> 
  filter(seqnames == 'chr4_h1', start > 27090375, end < 28361065)

# let's look at one gap in chr4 first
# get index of chr4 in subject
G_GR[seqnames(G_GR) == 'chr4_h1']

subjectHits(ol)

sat_GR[queryHits(ol)]

pairs <- findOverlapPairs(query = sat_GR, subject = G_GR,
                          ignore.strand = TRUE,
                          type = 'within') # only completely contained query objects are counted
ans <- pairs
mcols(ans)$overlap_width <- width(pintersect(ans, ignore.strand = TRUE))
ans

gap_TE$gapID <- gap_df$gapID[subjectHits(ol)]






library(HelloRanges)
bedtoo
# G_GR is the subject
subsetByOverlaps(G_GR, GR)

subsetByOverlaps(TE_S3h1_chr, gaps(gr_S3h2))
## TEs in the gaps
gap_TE <- subsetByOverlaps(TE_S3h1_chr, gaps(gr_S3h2))
# overlap for each gap
ol <- findOverlaps(query = TE_S3h1_chr, subject = gaps(gr_S3h2))
gap_TE$color <- element_color[match(gap_TE$simple_class, names(element_color))]
# use overlap info to assign the zoom names
gap_TE$gapID <- gap_df$gapID[subjectHits(ol)]









# data$y=runif(nrow(data), min = 0.8, max = 1.2)
# circos.genomicTrack(
#   data, stack = TRUE, track.index = 3,
#   panel.fun = function(region, value, ...) {
#     circos.genomicPoints(region, value, col = 'red', border = NA,
#                          lwd = 4, y = 1, pch =20)
#   })

# Create a dataframe for the legend
legend_df <- data.frame(TYPE = type_labels, color = color_palette)

# Create the ggplot2 legend
ggplot(legend_df, aes(x = 1, y = TYPE, fill = TYPE)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = color_palette, guide = guide_legend(title = "TYPE")) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm")
  )



