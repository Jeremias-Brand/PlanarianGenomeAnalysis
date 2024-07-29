library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)

txdb <- makeTxDbFromGFF("annotation/schMedS3_h1_ENCODE_hybrid_agat_only_hconf.gff3", format="gff3", dataSource=NA, organism=NA, taxonomyId=NA, chrominfo=NULL, miRBaseBuild=NA)

df <- read.table(file = 'conserved_elements/smed_base_conservation_annotated_peaks.tsv', sep = '\t',header = TRUE) |> 
  dplyr::select(chrom:summit, element_type, peak)

elements <- df |> 
  group_by(element_type) |> 
  group_split()

element_names <- levels(factor(df$element_type))
names(elements) <- element_names

################################################################################
## make peak files for each element type
for (e in element_names) {
  write_tsv(x = elements[[e]], file = paste0('supporting_information/', e, '.bed'), col_names = FALSE)
}



library(purrr)
annotated_peaks <- lapply(element_names, function(x) {
  annotatePeak(paste0('supporting_information/', x, '.bed'),
               TxDb = txdb, level = "transcript", assignGenomicAnnotation = TRUE, genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"),
               annoDb = NULL, addFlankGeneInfo = FALSE, flankDistance = 5000, sameStrand = FALSE, ignoreOverlap = FALSE, ignoreUpstream = FALSE,
               ignoreDownstream = FALSE, overlap = "TSS", verbose = TRUE)
}
)
names(annotated_peaks) <- element_names

df_peaks <- lapply(annotated_peaks, as.data.frame)

################################################################################

## summarizes and simplify the annotation

################################################################################

sum_peaks <- lapply(names(df_peaks), function(name) {
  df <- df_peaks[[name]]
  df |>
    mutate(simple_annotation = case_when(
      grepl('Promoter', annotation) ~ 'Promoter',
      grepl('Exon', annotation) ~ 'Exon',
      grepl('Intron', annotation) ~ 'Intron',
      grepl('Distal', annotation) ~ 'Distal',
      grepl('Downstream', annotation) ~ 'Downstream',
      TRUE ~ 'ERROR'
    )) |>
    mutate(simple_annotation = factor(simple_annotation, levels = c('Promoter', 'Exon', 'Intron', 'Distal', 'Downstream'))) |>
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
  mutate(percentage = N / total * 100)


pp <- ggplot(final_df, aes(x = element_class, y = percentage, fill = simple_annotation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.4) +
  geom_text(aes(label = total, y = 100, fill = NULL), vjust = -0.5, position = position_dodge(width = 0.9), 
            data = final_df %>% group_by(element_class) %>% summarise(total = sum(N)), check_overlap = TRUE) +
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#D55E00")) +
  theme_minimal() +
  labs(y = "Percentage", x = "", 
       fill = "") +
  theme(text = element_text(size = 10),
        legend.position = 'bottom')


pdf('fig/element_annotation.pdf', width = 4, height = 4)
pp
dev.off()

ggplot(final_df, aes(x = element_class, y = percentage, fill = simple_annotation)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_text(aes(label = N), vjust = -0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#D55E00")) +
  theme_minimal() +
  labs(y = "Percentage (%)", x = "Element Class", 
       fill = "Annotation Type", 
       title = "Percentage of Simple Annotations per Element Class")
      
      
ppie <- ggplot(final_df, aes(x = element_class, y = percentage, fill = simple_annotation)) +
  geom_pie  (stat = "identity", position = "stack", width = 0.4) +
  geom_text(aes(label = total, y = 100, fill = NULL), vjust = -0.5, position = position_dodge(width = 0.9), 
            data = final_df %>% group_by(element_class) %>% summarise(total = sum(N)), check_overlap = TRUE) +
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#D55E00")) +
  theme_minimal() +
  labs(y = "Percentage", x = "", 
       fill = "") +
  theme(text = element_text(size = 10),
        legend.position = 'bottom')

pp +
  coord_polar("y", start=0)

final_df |> 
  dplyr::select(simple_annotation, element_class, N, percentage) |> 
  pivot_wider(names_from = simple_annotation, values_from = c(N, percentage), id_cols = element_class) |> 
  write_xlsx('supporting_information/element_annotation.xlsx')

all_annotated <- do.call(rbind, df_peaks) 


all_annotated$peak = paste(all_annotated$seqnames, all_annotated$start, all_annotated$end, sep = '_')

ann2 <- all_annotated |> 
  dplyr::select(V12, annotation:distanceToTSS) |> 
  dplyr::rename(peak = V12) |> 
  dplyr::select(-c(geneChr, geneStart, geneEnd, geneLength, geneStrand))

peak_file <- read.table('conserved_elements/smed_base_conservation_annotated_peaks.tsv', header = TRUE)

left_join(peak_file, ann2, by = 'peak') |> 
  dplyr::select(peak:element_type, annotation:distanceToTSS, everything()) |> 
  write.table('supporting_information/supporting_selection/smed_all_annotations.tsv',
             sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
