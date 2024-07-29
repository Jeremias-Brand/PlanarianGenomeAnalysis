# look at genespace output
library(data.table)
library(tidyverse)
library(GENESPACE)
library(rtracklayer)
library(RColorBrewer)
library(ggbeeswarm)

################################################################################

## replotting with the generated object

################################################################################

load('genespace/results/gsParams.rda',
      verbose = TRUE)

new_gsParam <- gsParam
#gsParam$paths
#names(gsParam$paths)

#new_gsParam$paths <- lapply(gsParam$paths, FUN = function(x) {
#  gsub('./',x)})

# to plot with bp length we need to set useOrder = FALSE
psch <- plot_riparian(
  gsParam = new_gsParam, useOrder = FALSE,
  genomeIDs = rev(c(
    'schMedS3h1',
    'schMedS3h2',
    'schPol2',
    'schNov1',
    'schLug1'
    )),
  refGenome = 'schMedS3h1',
  refChrCols = brewer.pal(n=6, name  = 'Dark2'),
  chrFill = 'grey50',
  pdfFile = './genespace/genespace_S3.pdf') +
  theme_minimal()

################################################################################
## The matching blocks between genomes
sel = c('schMedS3h2', 'schPol2', 'schNov1', 'schLug1')
block_coords <- fread('./genespace/results/blkCoords.txt')

# calculate syntenic block size
blk_sizes <- block_coords |>
  filter(genome1 == 'schMedS3h1',
         genome2 %in% sel,
         grepl('[Cc]hr', chr1) & grepl('[Cc]hr', chr2)) |>
  mutate(block_size = abs(endBp1 - startBp1))

# calculate summary statistics for each genome
blk_sizes |> 
  group_by(genome2) |> 
  summarise(blocks = n(), 
            median = median(block_size),
            sd=sd(block_size),
            min=min(block_size),
            max=max(block_size),
            total_length=sum(block_size)) |> 
  write.table(file= 'genespace/schmidtea_synteny_block_size.tsv',
              quote = FALSE, sep = '\t', row.names = FALSE)

# 
pdf('genespace/genespace_synteny_block_size.pdf', width = 7)
blk_sizes |>
  ggplot(aes(y=block_size/1e6,  x= forcats::fct_relevel(genome2, sel))) +
  geom_beeswarm(color = 'grey40', dodge.width = 0) +
  geom_boxplot(alpha=0.6, fill = 'grey', width = 0.4) + 
  #facet_wrap(~chr1, ncol = 1, scales = 'free') +
  scale_color_brewer(palette = 'Dark2') +
  scale_y_log10(breaks = c(0.1,0.25,0.5,1,2.5,5,10,25,50,100,200)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size=20)) +
  ylab('synteny block size (Mb)') +
  scale_x_discrete(label=c(
    'S. med h2',
    'S. pol',
    'S. nov',
    'S. lug'
  ))
dev.off()

# 
pdf('genespace/genespace_synteny_gene_number.pdf', width = 7)
blk_sizes |>
  mutate(genes_in_block = nHits1) |>
  ggplot(aes(y=genes_in_block,  x= fct_relevel(genome2, sel))) +
  geom_beeswarm(color = 'grey40', dodge.width = 0) +
  geom_boxplot(alpha=0.6, fill = 'grey', width = 0.4) + 
  #facet_wrap(~chr1, ncol = 1, scales = 'free') +
  scale_color_brewer(palette = 'Dark2') +
  scale_y_log10(breaks = c(6,10,15,25,50,100,200, 500, 1000, 2000, 3000)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size=20)) +
  ylab('Number of genes in synteny block') +
  scale_x_discrete(label=c(
    'S. med h2',
    'S. pol',
    'S. nov',
    'S. lug'
  ))
dev.off()

