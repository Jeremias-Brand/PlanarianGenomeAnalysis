genome_size_from_fai <- function(fai) {
  fai <- fread(fai)
  return(sum(fai$V2))
}

rev_x = function(x, xrange = CELL_META$xlim) {
  xrange[2] - x + xrange[1]
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
