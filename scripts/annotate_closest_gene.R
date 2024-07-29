library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)
library(purrr)
# devtools::install_github("ebbertd/chisq.posthoc.test")
library(chisq.posthoc.test)
library(readxl)

AnnotationPriority = c("Promoter", "5UTR",  "Exon", "Intron", "3UTR", "Downstream", "Intergenic")
PromoterDefinition = c(-3000,0)
PlottingOrder = c('<1kb TSS', '1-2kb TSS', '2-3kb TSS', '5UTR',
 'Exon', 'Intron', '3UTR','Downstream <=300bp',  'Distal intergenic')

GFF="annotation/schMedS3_h1_ENCODE_hybrid_agat_only_hconf.gff3"
system(paste0('gunzip ', GFF, '.gz'))
txdb <- makeTxDbFromGFF(GFF, format="gff3", dataSource=NA, organism=NA, taxonomyId=NA, chrominfo=NULL, miRBaseBuild=NA)
system(paste0('gzip ', GFF))


df <- read.table(file = './conserved_elements/schMedS3h1_based_full_annotated_smed_peaks.tsv', sep = '\t',header = TRUE) 

elements <- df |> 
  group_by(element_type) |> 
  group_split()

element_names <- levels(factor(df$element_type))
names(elements) <- element_names

################################################################################
## make peak files for each element type
for (e in element_names) {
  write_tsv(x = elements[[e]], file = paste0('peak_annotation/', e, '.bed'), col_names = FALSE)
}

annotated_peaks <- lapply(element_names, function(x) {
  annotatePeak(paste0('peak_annotation/', x, '.bed'),
  tssRegion = PromoterDefinition,
  genomicAnnotationPriority = AnnotationPriority,
  TxDb = txdb, level = "transcript", assignGenomicAnnotation = TRUE, 
  annoDb = NULL, addFlankGeneInfo = FALSE, flankDistance = 5000, sameStrand = FALSE,
   ignoreOverlap = FALSE, ignoreUpstream = FALSE,
  ignoreDownstream = FALSE, overlap = "all", verbose = TRUE)
}
)
names(annotated_peaks) <- element_names

df_peaks <- lapply(annotated_peaks, as.data.frame)

OutputTable = bind_rows(df_peaks) |>
  mutate(simple_annotation = case_when(
        grepl('1kb', annotation) ~ '<1kb TSS',
        grepl('1-2kb', annotation) ~ '1-2kb TSS',
        grepl('2-3kb', annotation) ~ '2-3kb TSS',
        grepl('5\' UTR', annotation) ~ '5UTR',
        grepl('3\' UTR', annotation) ~ '3UTR',
        grepl('Exon', annotation) ~ 'Exon',
        grepl('Intron', annotation) ~ 'Intron',
        grepl('Distal', annotation) ~ 'Distal intergenic',
        grepl('Downstream', annotation) ~ 'Downstream <=300bp',
        TRUE ~ 'ERROR'
      )) |>
      rename(peak = V4,
      conservation_score = V5,
      element_type = V6,
      conservation_class = V7)
###########
# OUTPUT the annotation
OutputTable |> 
  write.table('./peak_annotation/schMedS3h1_by_class_closest.tsv', sep = '\t', row.names = FALSE)

################################################################################

## summarizes with more detail

################################################################################

sum_peaks <- lapply(names(df_peaks), function(name) {
  df <- df_peaks[[name]]
  df |>
    mutate(simple_annotation = case_when(
      grepl('1kb', annotation) ~ '<1kb TSS',
      grepl('1-2kb', annotation) ~ '1-2kb TSS',
      grepl('2-3kb', annotation) ~ '2-3kb TSS',
      grepl('5\' UTR', annotation) ~ '5UTR',
      grepl('3\' UTR', annotation) ~ '3UTR',
      grepl('Exon', annotation) ~ 'Exon',
      grepl('Intron', annotation) ~ 'Intron',
      grepl('Distal', annotation) ~ 'Distal intergenic',
      grepl('Downstream', annotation) ~ 'Downstream <=300bp',
      TRUE ~ 'ERROR'
    )) |>
    mutate(simple_annotation = factor(simple_annotation, levels = PlottingOrder)) |>
    group_by(simple_annotation) |>
    summarise(N = n()) |> 
    mutate(list_name = name) # Add the name of the list element as a column
  })

final_df <- bind_rows(sum_peaks) |> 
  mutate(element_class = factor(
    list_name, levels = c('putative_promoter', 'putative_enhancer', 'uncharacterized')))


# Calculate the total for each element_class
final_df <- final_df %>%
  group_by(element_class) %>%
  mutate(total = sum(N)) %>%
  ungroup()

# Calculate percentages
final_df <- final_df %>%
  mutate(percentage = N / total * 100) # |>
  # mutate(simple_annotation = factor(simple_annotation, levels = c('<1kb TSS', '1-2kb TSS', '2-3kb TSS', 'Exon', 'Intron','Downstream <=300bp',  'Distal intergenic', )))

col = rev(c('grey', '#C9B2D6', '#FE7F00', '#FCBF6F',  '#FB9A99', '#B2DE89', '#1F78B3', '#A5CEE3' ))
# okabe ito colors
col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
p2 <- ggplot(final_df, aes(x = element_class, y = percentage, fill = simple_annotation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = total, y = 100, fill = NULL), vjust = -0.5, position = position_dodge(width = 0.9), 
            data = final_df %>% group_by(element_class) %>% summarise(total = sum(N)), check_overlap = TRUE) +
  scale_x_discrete(labels = c('promoter', 'enhancer', 'unchar.')) +
  scale_fill_manual(values = col) +
  theme_minimal() +
  labs(y = "Percentage of peaks", x = "", 
       fill = "") +
  theme(text = element_text(size = 10),
        legend.position = 'right')


pdf('./peak_annotation/detailed_element_annotation.pdf', width = 4, height = 3)
p2
dev.off()

################################################################################

## STATISTICAL TEST if classes are enriched

################################################################################
m2 <- xtabs( ~ simple_annotation + element_type, data = OutputTable)
chisq.test(m2)
chisq.posthoc.test(m2, method = 'bonferroni', round = 1000) |> 
write.csv('./peak_annotation/schMedS3h1_distribution_by_element_type_chi_posthoc.csv')


################################################################################

## OTHER SPECIES

################################################################################

gffs = list(
  schMedS3h1 = "./annotation/schMedS3_h1_ENCODE_hybrid_agat_only_hconf.gff3",
  schPol2 = "./annotation/schPol2_ENCODE_hybrid_longest_coding.gff3",
  schLug1 = "./annotation/schLug1_ENCODE_hybrid_longest_coding.gff3",
  schNov1 = "./annotation/schNov1_ENCODE_hybrid_longest_coding.gff3"
)

annotate_them <- function(GFF, NAME){
  system(paste0('gunzip ', GFF, '.gz'))
  txdb <- makeTxDbFromGFF(GFF, format="gff3", dataSource=NA, organism=NA, taxonomyId=NA, chrominfo=NULL, miRBaseBuild=NA)
  system(paste0('gzip ', GFF))

  df <- read.table(file = paste0('./conserved_elements/ATAC_peak_files/', NAME, '_full_peak.tab'), sep = '\t',header = FALSE) 
  write_tsv(x = df, file = paste0('./conserved_elements/ATAC_peak_files/', NAME, '_full_peak.bed'), col_names = FALSE)

  print('annotating')
  ann <- annotatePeak(paste0('./conserved_elements/ATAC_peak_files/', NAME, '_full_peak.bed'),
  tssRegion = PromoterDefinition,
  genomicAnnotationPriority = AnnotationPriority,
  TxDb = txdb, level = "transcript", assignGenomicAnnotation = TRUE, 
  annoDb = NULL, addFlankGeneInfo = FALSE, flankDistance = 5000, sameStrand = FALSE,
  ignoreOverlap = FALSE, ignoreUpstream = FALSE,
  ignoreDownstream = FALSE, overlap = "all", verbose = TRUE)

  return(ann)
}

annotation_list <- map2(gffs, names(gffs), annotate_them)
schmidtea_peaks <- lapply(annotation_list, as.data.frame)


schmidtea_peaks2 <- lapply(names(schmidtea_peaks), function(name) {
  df <- schmidtea_peaks[[name]]
  df$species <- name
  df |>
  mutate(simple_annotation = case_when(
      grepl('1kb', annotation) ~ '<1kb TSS',
      grepl('1-2kb', annotation) ~ '1-2kb TSS',
      grepl('2-3kb', annotation) ~ '2-3kb TSS',
      grepl('5\' UTR', annotation) ~ '5UTR',
      grepl('3\' UTR', annotation) ~ '3UTR',
      grepl('Exon', annotation) ~ 'Exon',
      grepl('Intron', annotation) ~ 'Intron',
      grepl('Distal', annotation) ~ 'Distal intergenic',
      grepl('Downstream', annotation) ~ 'Downstream <=300bp',
      TRUE ~ 'ERROR'
    )) |>
    mutate(simple_annotation = factor(simple_annotation, levels = PlottingOrder))
  })


###

################################################################################

## OUTPUT

################################################################################
names(schmidtea_peaks2) <- names(schmidtea_peaks)

for (i in 1:length(schmidtea_peaks2)) {
  schmidtea_peaks2[[i]] |>
  dplyr::rename(peak = V4) |> 
   write.table( 
   paste0('./peak_annotation/', names(schmidtea_peaks2)[[i]], '_closest.tsv'),
    sep = '\t', row.names = FALSE)
}
################################################################################



names(schmidtea_peaks2) <- names(schmidtea_peaks)
sum_schmidtea_peaks <- lapply(names(schmidtea_peaks2), function(name) {
  df <- schmidtea_peaks2[[name]]
  df |>
    group_by(simple_annotation) |>
    summarise(N = n()) |> 
    mutate(list_name = name) # Add the name of the list element as a column

  })


sch_df  <- bind_rows(sum_schmidtea_peaks)

sch_df  <- sch_df  %>%
  group_by(list_name) %>%
  mutate(list_name = factor(list_name, levels = c(
    'schMedS3h1',
    'schPol2',
    'schLug1',
    'schNov1'))) |>
  mutate(total = sum(N)) %>%
  ungroup()

# Calculate percentages
sch_df <- sch_df %>%
  mutate(percentage = N / total * 100)

# col = rev(c('#C9B2D6', '#FE7F00', '#FCBF6F', '#FB9A99', '#B2DE89', '#1F78B3', '#A5CEE3'))
sch_p <- ggplot(sch_df, aes(x = list_name, y = percentage, fill = simple_annotation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = total, y = 100, fill = NULL), vjust = -0.5, position = position_dodge(width = 0.9), 
            data = sch_df  %>% group_by(list_name) %>% summarise(total = sum(N)), check_overlap = TRUE) +
  scale_fill_manual(values = col) +
  theme_minimal() +
  labs(y = "Percentage of peaks", x = "", 
       fill = "") +
  theme(text = element_text(size = 10),
        legend.position = 'right')



pdf('./peak_annotation/schmidtea_element_annotation.pdf', width = 4, height = 3)
sch_p
dev.off()


################################################################################
tbl = bind_rows(schmidtea_peaks2)
m2 <- xtabs( ~ simple_annotation + species, data = tbl)
sum(m2[])
m3 <- m2[rowSums(m2) != 0,]
chisq.test(m3)
chisq.posthoc.test(m3, method = 'bonferroni', round = 1000) |> 
write.csv('./peak_annotation/Schmidtea_annotation_by_species_chi_posthoc.csv')




#############

# Annotation chip peaks
file  <- './supporting_information/Additional File 2.xlsx'
atac  <- read_xlsx(file, sheet = "TabS1_schMedS3h1_ATAC")
me  <- read_xlsx(file, sheet = "TabS2_schMedS3h1_wtH3K4me3")
ac  <- read_xlsx(file, sheet = "TabS3_schMedS3h1_wtH3K27ac")


peak_files <- list(atac, me, ac)

################################################################################
## make peak files for each element type


annotated_peaks <- lapply(peak_files, function(x) {
    annotatePeak(makeGRangesFromDataFrame(x),
    tssRegion = PromoterDefinition,
    genomicAnnotationPriority = AnnotationPriority,
    TxDb = txdb, level = "transcript", assignGenomicAnnotation = TRUE, 
    annoDb = NULL, addFlankGeneInfo = FALSE, flankDistance = 5000, sameStrand = FALSE,
    ignoreOverlap = FALSE, ignoreUpstream = FALSE,
    ignoreDownstream = FALSE, overlap = "all", verbose = TRUE)
  }
)

df_types <- lapply(annotated_peaks, as.data.frame)



selected_columns <- lapply(df_types, function(df) {
    df[, c("geneId", "distanceToTSS")]
})

selected_columns[[1]]$mark  <- 'ATACseq'
selected_columns[[2]]$mark  <- 'H3K4me3'
selected_columns[[3]]$mark  <- 'H3K27ac'

plotting_df <- do.call(rbind, selected_columns)

breaks <- seq(-5e5, 5e5, by = 1e5)  # Modify as needed
labels <- paste0(breaks / 1e3, " kb")  # Converting units to kb for readability

pdf('element_distance_TSS_smed.pdf', width = 3.5, height = 7)
ggplot(plotting_df, aes(x = distanceToTSS)) +
geom_histogram(aes(fill = mark), color = 'black', binwidth = 50000) +
scale_y_log10() +
#xlim(c(-5e5,5e5)) +
scale_x_continuous(breaks = breaks, labels = labels) +
facet_wrap(~mark, nrow = 3) +
scale_fill_manual(values = c("#E69F00",  "#0072B2", "#CC79A7")) +
theme_minimal() +
theme(axis.title = element_blank(),
axis.text = element_text(size = 10, angle = 90),
facet.text = element_text(size = 8),
legend.position = 'none')

dev.off()

plotting_df |> 
group_by(mark) |> 
summarise(
  median = median(distanceToTSS),
  SD = sd(distanceToTSS)
)

# #
# names(annotated_peaks) <- element_names


# col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")


# df <- df_types[[2]]
# break_size = 5e4

# brs <- c(-Inf, seq(-1e6, 1e6, break_size), Inf)
# # adjust breaks with prior knowledge for vizualization
# brs <- c(-Inf, seq(-350e5, 250e5, break_size), Inf)

# lb <- paste0(seq(-350e5, 250e5 - break_size, break_size) / 1e3, 'kb')
# lb2 <- ifelse(grepl('0', lb), lb, '')
# lb2 = lb
# lab  <- c('<-1Mb', lb2, '>1Mb')

# p3 <- df |> 
#   mutate(distTSSbins = cut(distanceToTSS, breaks = brs, labels = lab)) |>
#   group_by(distTSSbins) |>
#   summarise(N = n()) |>
#   ggplot(aes(x = distTSSbins, y = N)) +
#   geom_bar(stat = 'identity') +
#   scale_y_log10() +
#   #scale_x_discrete(labels = ifelse(grepl('0', lb), lb, '')) +
#   theme_classic() +
#   theme(axis.text = element_text(angle = 90, size = 20)) 

# summary(df$distanceToTSS)


# brs <- c(-Inf, -1e6, -1e5, -1e4, -1e3, -100, -10, 0, 10, 100, 1e3, 1e4, 1e5, 1e6, Inf)





# ggplot(plotting_df, aes(x = mark, y = distanceToTSS)) +
# geom_violin(aes(fill = mark), color = 'black', binwidth = 50000) +
# scale_y_log10() +
# #xlim(c(-5e5,5e5)) +
# scale_fill_manual(values = c("#E69F00", "#CC79A7", "#0072B2")) +
# theme_classic() +
# theme(axis.title = element_blank(),
# axis.text = element_text(size = 10))




# bin_width = 50000

# ggplot(plotting_df, aes(x = distanceToTSS)) +
# geom_histogram(aes(fill = mark), color = 'black', binwidth = bin_width) +
# scale_y_log10()  +
#   geom_text(stat = 'bin', aes(y = ..count.., label = ..count..), 
#             binwidth = bin_width, vjust = -0.5, hjust = -0.1, 
#             position = position_nudge(x = -bin_width / 2))






# lapply(df_types, '[', c("geneId", "distanceToTSS"))

# library(cowplot)
# plots = list()
# col = c("#E69F00", "#CC79A7", "#0072B2")
# for (i in 1:3){
# df <- df_types[[i]]
# plots[[i]]  <- ggplot(df, aes(x = distanceToTSS)) +
# geom_histogram(fill = col[i], color = 'black', binwidth = 50000) +
# scale_y_log10() +
# xlim(c(-5e5,5e5)) +
# theme_classic() +
# theme(axis.title = element_blank(),
# axis.text = element_text(size = 10))
# print(summary(df$distanceToTSS))
# }


# plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 1, align = 'hv', axis = 'tblr')
