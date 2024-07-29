# mamba install -y r-ggalluvial r-ggridges 

library(data.table)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(forcats)
library(cowplot)
library(stringr)
library(ggridges)
library(RColorBrewer)
library(ggplot2)

################################################################################

## look at the output

################################################################################

smed_based_result <- readRDS('./schMedS3h1_based_full_result.rds')


dfs <- lapply(smed_based_result, '[[', 1)
beds <- lapply(smed_based_result, '[[', 2)


source_bed = './ATAC_peak_files/schMedS3h1_full_peak.tab'
df_smed_peaks <- read.table(source_bed)
names(df_smed_peaks) <- c('chrom', 'start', 'end', 'name')

################################################################################

## let's deal with the df later

################################################################################

big_bed <- do.call(rbind, beds)
big_bed$species <- gsub('.[0-9]+$', '', row.names(big_bed))

new_df <- big_bed %>%
  mutate(peak_type = case_when(
    summit_OL_count == 1 & atac_OL_count == 1 ~ "peak_strict",
    (summit_OL_count == 1 & atac_OL_count >= 1) | 
      (summit_OL_count >= 1 & atac_OL_count == 1) | 
      (summit_OL_count > 1 & atac_OL_count > 1) ~ "peak_strict_multiOL",
    summit_OL_count == 0 & atac_OL_count == 1 ~ "peak_nS",
    summit_OL_count == 0 & atac_OL_count > 1 ~ "peak_nS_multiOL",
    summit_OL_count == 1 & atac_OL_count == 0 ~ "region",
    summit_OL_count > 1 & atac_OL_count == 0 ~ "region_multiOL",
    summit_OL_count == 0 & atac_OL_count == 0 ~ "region_nS",
    TRUE ~ "unknown"
  )) |> 
  mutate(simple_peak_type = case_when(
    peak_type %in% c("peak_strict", "peak_strict_multiOL") ~ 'cons_peak',
    peak_type %in% c("peak_nS", "peak_nS_multiOL") ~ 'peak_nS',
    peak_type %in% c("region_nS") ~ 'region_nS',
    TRUE ~ 'cons_region',
  ))


new_df |> 
  group_by(simple_peak_type, species) |> 
  summarise(N = n())

new_df <- new_df |> 
  mutate(simple_peak_type = factor(simple_peak_type,
                                   levels = c("cons_peak", "peak_nS", "cons_region",  "region_nS")))
                                   

                       
# we now need to summarise these data per query peak
df <- new_df |> 
  group_by(name.distinct, species, simple_peak_type) |> 
  summarise(N = n()) |> 
  filter(simple_peak_type %in% c('cons_peak', 'cons_region'))

df2 <- df %>%
  group_by(name.distinct, species) %>%
  summarise(simple_peak_type = case_when(
    any(simple_peak_type == "cons_peak") ~ "cons_peak",
    TRUE ~ 'cons_region'
  )) 


bi_tmp <- df2 
bi_tmp$binary_status <-case_when(
  df2$simple_peak_type == 'cons_peak' ~ 10,
  df2$simple_peak_type == 'cons_region' ~ 01,
  TRUE ~ 00
)

################################################################################
## add the statuses creating a binary indicator of conservation
bi <- bi_tmp |>
  group_by(name.distinct) |>
  summarise(binary_sum = sum(binary_status))

bi <- left_join(df_smed_peaks, bi, by = c('name' = 'name.distinct')) |> 
  mutate(binary_sum = ifelse(is.na(binary_sum), 00, binary_sum))

write.table(bi, file = './schMedS3h1_base_full_conservation_scoring_binary.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

# bi <- read.table('./conserved_elements/schMedS3h1_base_full_conservation_scoring_binary.tsv', header = TRUE)
bi |>
group_by(binary_sum) |>
  summarise(N = n()) |> 
  mutate(pcnt = N/sum(N) * 100)


################################################################################

## add the enhancer and conditions data

################################################################################

# enhancer: uniq_chr_Smed_combined_peaks_wt_H3K27ac_noH3K4me3_summits_400bp.tab
# promoter: uniq_chr_Smed_combined_peaks_wt_H3K27ac_withH3K4me3_summits_400bp.tab
# all:      chr_Smed_combined_summits_400bp.tab

promoter <- read.table('./../final_putative_enhancers/uniq_chr_Smed_combined_peaks_wt_H3K27ac_withH3K4me3_summits_400bp.tab')
enhancer <- read.table('./../final_putative_enhancers/uniq_chr_Smed_combined_peaks_wt_H3K27ac_noH3K4me3_summits_400bp.tab')

bi$putative_promoter <- ifelse(bi$name %in% promoter$V4, 'YES', 'NO')
bi$putative_enhancer <- ifelse(bi$name %in% enhancer$V4, 'YES', 'NO')

bi  |> 
filter(putative_enhancer == 'YES' & putative_promoter == 'YES')  |> 
View()

qval = 0.05
# df <- read.csv('doc/supporting_information/BEDlike_files_comparative_genomics/Orthologous_accessible_regions/comprehensive_accessible_chromatin_all.csv')

################################################################################
## ATAC peaks classified as both putative enhancer (no K4me) and promoter have both signals and should be classified as promoters
df <- bi |> 
  mutate(
    element_type = case_when(
      putative_enhancer == 'YES' & putative_promoter == 'YES' ~ 'putative_promoter',
      putative_promoter == 'YES' ~ 'putative_promoter',
      putative_enhancer == 'YES'  ~ 'putative_enhancer',
      TRUE ~ 'uncharacterized'
    )
  ) |> 
  dplyr::rename(conservation_score = binary_sum) |> 
  mutate(conservation_class = case_when(
    conservation_score == 30 ~ 'highly_conserved',
    conservation_score != 0 ~ 'partially_conserved',
    TRUE ~ 'not_conserved'
    )) |>
  rename(peak = name) |>
  dplyr::select(-c(putative_enhancer, putative_promoter))


################################################################################

## plotting function

################################################################################
print('define plotting functions')
# Define a custom function to create the plot
create_plot <- function(data, group_var, cols, lab) {
  # Summarize data and calculate percentages
  summary_df <- data %>%
    filter(!is.na(conservation_score)) %>%
    group_by(!!sym(group_var), conservation_score) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / sum(count) * 100, .groups = 'drop')
  
  # Calculate total counts for each group_var
  total_counts <- summary_df %>%
    group_by(!!sym(group_var)) %>%
    summarise(total = sum(count), .groups = 'drop')
  
  # Create the plot
  plot <- summary_df %>%
    ggplot(aes(x = !!sym(group_var), y = percentage, fill = factor(conservation_score))) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(
      values = cols,
      name = "conservation type",
      labels = lab) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          axis.title = element_blank()) +
    geom_text(aes(label = paste0(round(percentage, 1))),
              position = position_stack(vjust = 0.3), size = 4, color = 'white') + 
    geom_text(data = total_counts, aes(x = !!sym(group_var), y = 100, label = total, fill = NULL),
              vjust = -0.5, size = 4, color = 'black', fill = NULL, check_overlap = TRUE)
  
  return(plot)
}
plot('using function')
p_smed_cons <- df |>
  ggplot(aes(x = 'all peaks', fill = factor(conservation_score))) +
  geom_bar(width = 0.8) +
  #geom_hline(yintercept = nrow(res_df)) +
  scale_fill_manual(
    values =
      colors <- c("grey60", "#fe8181", "#fe5757", "#cb2424", "#43c4ef", "#00b4ff", "#009eff", "#0066ff", "#004fff", "#000080"),
    name = "conservation type",
    labels = c("not conserved",
               "1 region",
               "2 region",
               "3 region",
               "1 peak",
               "1 peak, 1 region",
               "1 peak, 2 region",
               "2 peak",
               "2 peak, 1 region",
               "3 peak")) +
  theme_minimal() +
  xlab('') +
  ylab('') +
  # Add percentages to the plot
  geom_text(stat = "count", aes(label = paste0(round((..count..)/sum(..count..)*100, 1))),
            position = position_stack(vjust = 0.5), size = 4, color = 'white') +
  geom_text(y = 55585, label = '55585',
              vjust = -0.5, size = 4, color = 'black', fill = NULL, check_overlap = TRUE)

l1 <- get_legend(p_smed_cons + theme(legend.position = 'right'))

# 

print('using function')
p_highly_cons <- df |>
  filter(conservation_score == 30) |> 
  mutate(element_type = factor(element_type, levels = c(
    'putative_promoter', 'putative_enhancer', 'uncharacterized'
  ))) |>
  ggplot(aes(x = 'highly_conserved_peaks', fill = factor(element_type))) +
  geom_bar(width = 0.8) +
  #geom_hline(yintercept = nrow(res_df)) +
  scale_fill_manual(
    values =
      colors <- c("#CC79A7", "#0072B2", "#E69F00"),
    name = "conservation type",
    labels = c('promoter', 'enhancer', 'uncharacterized')) +
  theme_minimal() +
  xlab('') +
  ylab('') +
  # Add percentages to the plot
  geom_text(stat = "count", aes(label = paste0(round((..count..)/sum(..count..)*100, 1))),
            position = position_stack(vjust = 0.5), size = 4, color = 'white') +
  geom_text(y = 7567, label = '7567',
              vjust = -0.5, size = 4, color = 'black', fill = NULL, check_overlap = TRUE)

l2 <- get_legend(p_highly_cons + theme(legend.position = 'right'))
##

print('using function')
# Define colors and labels
cols <- c("grey60", "#fe8181", "#fe5757", "#cb2424", "#43c4ef", "#00b4ff", "#009eff", "#0066ff", "#004fff", "#000080")
lab <- c("not conserved", "1 region", "2 region", "3 region", "1 peak", "1 peak, 1 region", "1 peak, 2 region", "2 peak", "2 peak, 1 region", "3 peak")

# Create plots for different grouping variables
element_type_plot <- df |> 
  mutate(element_type = factor(element_type, levels = c('putative_promoter', 'putative_enhancer', 'uncharacterized'))) |> 
           create_plot("element_type", cols, lab)
         
# element_type_plot 


p1 <- plot_grid(
          p_smed_cons + theme(legend.position = 'none'),
          element_type_plot + theme(legend.position = 'none'), 
          labels = c('A', 'B'),
          rel_widths = c(1, 3), 
          align = 'hv')


p2 <- plot_grid(p1, l1, 
                rel_widths = c(1,0.2))

# Show the plots
pdf('./conservation_by_annotation.pdf')
print(p2)
dev.off()


# Plot with percentage of conserved that has a ChIP mark
pdf('./conservation_by_annotation_with_element_type.pdf')
p1 <- plot_grid(
          p_smed_cons + theme(legend.position = 'none'),
          p_highly_cons+ theme(legend.position = 'none'), 
          element_type_plot + theme(legend.position = 'none'), 
          labels = c('A', 'B', 'C'),
          rel_widths = c(1, 1, 3), 
          align = 'hv', nrow = 1 )


p2 <- plot_grid(l1, l2, nrow = 2, 
                rel_heights = c(2,1))
p3 <- plot_grid(p1, p2, 
                rel_widths = c(1,0.2))
print(p3)
dev.off()
################################################################################
## output
write.table(df, file = './schMedS3h1_based_full_annotated_smed_peaks.tsv',
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################

## output bed files for each liftover species

################################################################################
################################################################################
## For the bed files we would like to enhance it with the categorization

beds <- lapply(smed_based_result, '[[', 2)
filtered_bi <- df |> 
  dplyr::select(peak, conservation_class, conservation_score) 

joined_beds <- lapply(beds, function(x) {
  left_join(x, filtered_bi, by = c('name.distinct' = 'peak')) |> 
    dplyr::select(-grouping) |> 
    dplyr::rename(smed_peak = name.distinct)
})

# Function to write each element of the list to a file
write_list_elements <- function(element_name, element_data) {
  filename <- paste0('./smed_base_full_conservation_', element_name, "_coordinates.tsv") # Set the filename using the element name
  write.table(element_data, file = filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# write to output
for (i in 1:length(joined_beds)) {
  write_list_elements(names(joined_beds)[[i]], joined_beds[[i]])
}



cols = c('#4281A6', '#1034A6',  '#007FFF',
         '#FF033E', '#F9B5D0', 'grey60', 'black')


pdf('fig/peak_conservation_supplement.pdf')
ggplot(new_df, aes(x = width, y = simple_peak_type, fill = simple_peak_type)) + 
  geom_boxplot(
    position = position_nudge(y=-.1),
    width = .12, 
    outlier.color = NA 
  ) +
  geom_density_ridges(alpha = 0.4, scale = 0.8, size = 0.25) +
  facet_wrap(~ species, scales = "free_y", ncol = 1) +
  labs(x = "Size of liftover",
       y = "",
       fill = "Peak Type") +
  scale_fill_manual(values = cols[2:5]) +
  scale_y_discrete(labels = c('cons. peak', 'peak no summit', 'cons. region',  'region no summit')) +
  theme_minimal() +
  theme(legend.position = 'none') 
dev.off()

p_regions <- ggplot(new_df, aes(x = width, y = simple_peak_type, fill = simple_peak_type)) + 
  geom_boxplot(
    position = position_nudge(y=-.1),
    width = .12, 
    outlier.color = NA 
  ) +
  geom_density_ridges(alpha = 0.4, scale = 0.8, size = 0.25) +
  # geom_vline(xintercept = 400) +
  labs(
    x = "Size of liftover",
    y = "",
    fill = "Peak Type") +
  scale_fill_manual(values = cols[2:5]) +
  scale_y_discrete(labels = c('cons. peak', 'peak no summit', 'cons. region',  'region no summit')) +
  scale_x_continuous(breaks = seq(0,1000, 200)) +
  theme_minimal() +
  theme(legend.position = 'none') 
pdf('fig/peak_conservation.pdf')

p_regions
dev.off()

################################################################################
## summary for a supplemental plot
new_df |> 
  group_by(species, simple_peak_type) |> 
  summarise(
    Median = median(width, na.rm = TRUE),
    Std_Dev = sd(width, na.rm = TRUE),
    Min = min(width, na.rm = TRUE),
    Max = max(width, na.rm = TRUE)
  ) |> 
  write.csv('./liftover_size_by_species_and_peak_type.tsv')

new_df |> 
  group_by(simple_peak_type) |> 
  summarise(
    Median = median(width, na.rm = TRUE),
    Std_Dev = sd(width, na.rm = TRUE),
    Min = min(width, na.rm = TRUE),
    Max = max(width, na.rm = TRUE)
  ) |> 
  write.csv('./liftover_size_by_peak_type.tsv')

################################################################################

## plotting

################################################################################
# 
devtools::install_github("ebbertd/chisq.posthoc.test")
## overall plot
#install.packages("devtools")
#devtools::install_github("ebbertd/chisq.posthoc.test")
library(chisq.posthoc.test)

tbl2 <- xtabs( ~ conservation_class + element_type, data = df)
round(tbl2/43900,2)

# over all
round(100*xtabs( ~ conservation_class, data = df)/43900,1)
xtabs( ~ conservation_class + element_type, data = df)


m2 <- xtabs( ~ conservation_class + element_type, data = df)
chisq.test(m2)
chisq.posthoc.test(m2, method = 'bonferroni', round = 4) |> 
write.csv('./chi_posthoc.csv')
