library(data.table)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(forcats)
library(cowplot)
library(stringr)
library(ggridges)
library(RColorBrewer)

################################################################################

## look at the output

################################################################################

smed_based_result <- readRDS('conserved_elements/smed_based_result.rds')


dfs <- lapply(smed_based_result, '[[', 1)
beds <- lapply(smed_based_result, '[[', 2)


source_bed = 'WGA/Smed_based/Smed_CRE_IDR/summit/modif_chr_SM_SMX_SMT_SMH_summits_400bp.tab'
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


write.table(bi, file = 'conserved_elements/smed_base_conservation_scoring_binary.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################

## add the enhancer and conditions data

################################################################################
qval = 0.05
df <- read.csv('doc/supporting_information/BEDlike_files_comparative_genomics/Orthologous_accessible_regions/comprehensive_accessible_chromatin_all.csv')

################################################################################
## ATAC peaks classified as both putative enhancer (no K4me) and promoter have both signals and should be classified as promoters
exp_df <- df |> 
  mutate(
    element_type = case_when(
      putative_enhancer == 'YES' & proximal_CRE == 'YES' ~ 'putative_promoter',
      proximal_CRE == 'YES' ~ 'putative_promoter',
      putative_enhancer == 'YES'  ~ 'putative_enhancer',
      TRUE ~ 'uncharacterized'
    )
  ) |> 
  dplyr::rename(summit = peak,
                peak = X) |> 
  # remove functional data information.
  dplyr::select(-c(putative_enhancer, proximal_CRE), -contains('p.value'), -contains('_Fold'), -contains('Conc_'), -contains('_FDR'))

################################################################################
## output
write.table(exp_df, file = 'final_putative_enhancers/annotated_smed_peaks.tsv',
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)


################################################################################
## Merge with conservation data
bi_exp_df <- left_join(exp_df, dplyr::select(bi, name, binary_sum),  by = c('peak' = 'name')) |> 
  dplyr::rename(conservation_score = binary_sum) |> 
  mutate(conservation_class = case_when(
    conservation_score == 30 ~ 'highly_conserved',
    conservation_score != 0 ~ 'partially_conserved',
    TRUE ~ 'not_conserved'
    )) |> 
  dplyr::select(peak:summit, element_type, conservation_class, conservation_score, everything())
  
bi_exp_df[bi_exp_df == ""] <- NA
################################################################################
## output
write.table(bi_exp_df, file = 'conserved_elements/smed_base_conservation_annotated_peaks.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################

## output bed files for each liftover species

################################################################################

################################################################################
## for the classification we can combine all output
res_lst <- list()

types <- lapply(smed_based_result, '[[', 1)
type_df <- do.call(rbind, types)

write.table(type_df, 'conserved_elements/smed_base_conservation_type_by_species.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################
## For the bed files we would like to enhance it with the categorization

beds <- lapply(smed_based_result, '[[', 2)
filtered_bi <- bi_exp_df |> 
  dplyr::select(peak, conservation_class, conservation_score) 

joined_beds <- lapply(beds, function(x) {
  left_join(x, filtered_bi, by = c('name.distinct' = 'peak')) |> 
    dplyr::select(-grouping) |> 
    dplyr::rename(smed_peak = name.distinct)
})

# Function to write each element of the list to a file
write_list_elements <- function(element_name, element_data) {
  filename <- paste0('conserved_elements/smed_base_conservation_', element_name, "_coordinates.tsv") # Set the filename using the element name
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
  write.csv('conserved_elements/liftover_size_by_species_and_peak_type.tsv')

new_df |> 
  group_by(simple_peak_type) |> 
  summarise(
    Median = median(width, na.rm = TRUE),
    Std_Dev = sd(width, na.rm = TRUE),
    Min = min(width, na.rm = TRUE),
    Max = max(width, na.rm = TRUE)
  ) |> 
  write.csv('conserved_elements/liftover_size_by_peak_type.tsv')

################################################################################

## plotting

################################################################################
## overall plot

p_smed_cons <- bi |>
  ggplot(aes(x = 'all peaks', fill = factor(binary_sum))) +
  geom_bar(width = 0.6) +
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
  geom_text(stat = "count", aes(label = paste0(round((..count..)/sum(..count..)*100), "%")),
            position = position_stack(vjust = 0.5), size = 4, color = 'white') 

l1 <- get_legend(p_smed_cons + theme(legend.position = 'right'))

################################################################################

## plotting function

################################################################################

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
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              position = position_stack(vjust = 0.3), size = 4, color = 'white') + 
    geom_text(data = total_counts, aes(x = !!sym(group_var), y = 100, label = total, fill = NULL),
              vjust = -0.5, size = 4, color = 'black', fill = NULL, check_overlap = TRUE)
  
  return(plot)
}

# Define colors and labels
cols <- c("grey60", "#fe8181", "#fe5757", "#cb2424", "#43c4ef", "#00b4ff", "#009eff", "#0066ff", "#004fff", "#000080")
lab <- c("not conserved", "1 region", "2 region", "3 region", "1 peak", "1 peak, 1 region", "1 peak, 2 region", "2 peak", "2 peak, 1 region", "3 peak")

# Create plots for different grouping variables
element_type_plot <- bi_exp_df |> 
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
# p2

# proximal_CRE_plot <- create_plot(bi_exp_df, "proximal_CRE", cols, lab)
# ht_group_plot <- create_plot(bi_exp_df, "ht_group", cols, lab)
# xray_group_plot <- create_plot(bi_exp_df, "xray_group", cols, lab)

# Show the plots
pdf('fig/conservation_by_annotation.pdf')
print(p2)
# print(proximal_CRE_plot)
# print(ht_group_plot)
# print(xray_group_plot)
dev.off()


#install.packages("devtools")
#devtools::install_github("ebbertd/chisq.posthoc.test")
library(chisq.posthoc.test)

tbl2 <- xtabs( ~ conservation_class + element_type, data = bi_exp_df)
round(tbl2/43900,2)


tidy(tbl2)

# over all
round(100*xtabs( ~ conservation_class, data = bi_exp_df)/43900,1)
xtabs( ~ conservation_class + element_type, data = bi_exp_df)


m2 <- xtabs( ~ conservation_class + element_type, data = bi_exp_df)
chisq.test(m2)
chisq.posthoc.test(m2, method = 'bonferroni', round = 4) |> 
write_csv('conserved_elements/chi_posthoc.csv')
################################################################################

## compound figure

################################################################################

create_plot2 <- function(data, group_var, cols, lab) {

  # Summarize data and calculate percentages
  summary_df <- data %>%
    group_by(!!sym(group_var), binary_sum) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / sum(count) * 100)
  
  # Create the plot
  plot <- summary_df %>%
    ggplot(aes(x = !!sym(group_var), y = percentage, fill = factor(binary_sum))) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(
      values = cols,
      name = "conservation type",
      labels = lab) +
    ylab('') +
    scale_y_continuous(expand = c(0, 4)) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x  = element_text(size = 14, color = 'black', angle = 45, hjust = 1.1),
          axis.text.y  = element_text(size = 14, color = 'black')) 
  
  return(plot)
}

putative_enhancer_plot <- create_plot2(bi_exp_df, "putative_enhancer", cols, lab)
proximal_CRE_plot <- create_plot2(bi_exp_df, "proximal_CRE", cols, lab)
ht_group_plot <- create_plot2(bi_exp_df, "ht_group", cols, lab)
xray_group_plot <- create_plot2(bi_exp_df, "xray_group", cols, lab)


l1 <- get_legend(putative_enhancer_plot + theme(legend.position = 'bottom', legend.title = element_blank()))
p1 <- plot_grid(putative_enhancer_plot + 
                  scale_x_discrete(labels = c('Other', 'Enhancer')) +
                  theme(legend.position = 'none', axis.title.x = element_blank()) , 
          proximal_CRE_plot  + 
            scale_x_discrete(labels = c('Other', 'CRE')) +
            theme(legend.position = 'none', axis.title.x = element_blank()),
          ht_group_plot  + 
            scale_x_discrete(labels = c('Head', 'Tail', 'Other')) +
            theme(legend.position = 'none', axis.title.x = element_blank()),
          xray_group_plot  + 
            scale_x_discrete(labels = c('WT', 'X-ray', 'Other')) +
            theme(legend.position = 'none', axis.title.x = element_blank()),
          nrow = 1 ,
          labels = LETTERS[1:4],
          align = 'h', axis = 'tb')

pdf('fig/conservation_by_annotation_main.pdf')
plot_grid(p1, l1, ncol = 1, rel_heights = c(1, 0.2))
dev.off()



################################################################################

## We would like to have a plot for Other, Enhancer, CRE

################################################################################


create_plot3 <- function(data, group_var, cols, lab, filt = 'YES', xlab) {
  # Summarize data and calculate percentages
  summary_df <- data %>%
    filter(!!sym(group_var) %in% filt) |> 
    group_by(!!sym(group_var), binary_sum) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / sum(count) * 100)
  
  # Create the plot
  plot <- summary_df %>%
    ggplot(aes(x = !!sym(group_var), y = percentage, fill = factor(binary_sum))) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(
      values = cols,
      name = "conservation type",
      labels = lab) +
    ylab('') +
    xlab(xlab) +
    scale_y_continuous(expand = c(0.05, 4), limits = c(0,100)) +
    theme_minimal() +
    theme(legend.position = 'none',
      panel.grid.major.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.text.y  = element_text(size = 14, color = 'black'))   +
    geom_text(aes(label = round(percentage,1)), position = position_stack(vjust = 0.5), size = 4, color = 'white') +
    geom_text(x=1 , y = 105, label = sum(summary_df$count), size = 5, check_overlap = TRUE)

  return(plot)
}

bi_exp_df$all = 'YES'

bi_exp_df$head <- ifelse(bi_exp_df$ht_group == "upT" | is.na(bi_exp_df$ht_group), 'NO', 'YES')
bi_exp_df$tail <- ifelse(bi_exp_df$ht_group == "upH" | is.na(bi_exp_df$ht_group), 'NO', 'YES')
bi_exp_df$wt <-   ifelse(bi_exp_df$xray_group == "upXray" | is.na(bi_exp_df$xray_group), 'NO', 'YES')
bi_exp_df$xray <-   ifelse(bi_exp_df$xray_group == "upWt" | is.na(bi_exp_df$xray_group), 'NO', 'YES')



pdf('fig/conservation_by_putative_element.pdf', height = 7, width = 4)
plot_grid(
create_plot3(bi_exp_df, "all", cols, lab, filt = 'YES', 'All') + theme(axis.text.y = element_blank(),
                                                                       panel.grid.major = element_blank() ,
                                                                       panel.grid.minor = element_blank()),
create_plot3(bi_exp_df, "putative_enhancer", cols, lab, filt = 'YES', 'putative enhancer') + theme(axis.text.y = element_blank(),
                                                                                                   panel.grid.major = element_blank() ,
                                                                                                   panel.grid.minor = element_blank()),
create_plot3(bi_exp_df, "proximal_CRE", cols, lab, filt = 'YES', 'proximal CRE') + theme(axis.text.y = element_blank(),
                                                                                         panel.grid.major = element_blank() ,
                                                                                         panel.grid.minor = element_blank()),
nrow = 1
)
dev.off()


pdf('fig/conservation_by_functional.pdf', height = 7, width = 6.66)
plot_grid(
  create_plot3(bi_exp_df, "all", cols, lab, filt = 'YES', 'All') + theme(axis.text.y = element_blank(),
                                                                         panel.grid.major = element_blank() ,
                                                                         panel.grid.minor = element_blank()),
  create_plot3(bi_exp_df, "head", cols, lab, filt = 'YES', 'Head') + theme(axis.text.y = element_blank(),
                                                                           panel.grid.major = element_blank() ,
                                                                           panel.grid.minor = element_blank()),
  create_plot3(bi_exp_df, "tail", cols, lab, filt = 'YES', 'Tail') + theme(axis.text.y = element_blank(),
                                                                           panel.grid.major = element_blank() ,
                                                                           panel.grid.minor = element_blank()),
  create_plot3(bi_exp_df, "wt", cols, lab, filt = 'YES', 'WT') + theme(axis.text.y = element_blank(),
                                                                       panel.grid.major = element_blank() ,
                                                                       panel.grid.minor = element_blank()),
  create_plot3(bi_exp_df, "xray", cols, lab, filt = 'YES', 'Xray') + theme(axis.text.y = element_blank(),
                                                                           panel.grid.major = element_blank() ,
                                                                           panel.grid.minor = element_blank()),
  nrow = 1
)
dev.off()


################################################################################

## To get a statistical we apply a bootstrapping approach

################################################################################

# use permutation test
permutation_test <- function(df, var, n_permutations = 10000, outfile) {
  ## DEBUG
  df = testing_df
  var = 'all_v_uncharacterized'
  
  # we are only testing for the 3 conserved peaks group
  test_df <-  df %>%
    as.data.frame() %>%
    mutate(cons = ifelse(binary_sum == 30, 1, 0)) 
  # Calculate the actual difference in the fraction of conserved between the putative enhancer groups
  all_frac <- test_df[,var]
  actual_df <- test_df |> 
    group_by(!!sym(var)) |> 
    summarise(N= n(),
              cons_sum = sum(cons),
              frac= cons_sum / n()) 
  # get the actual difference our test statistic
  actual_diff <-   actual_df |> 
    summarise(delta_frac = diff(frac)) |> 
    pull(delta_frac)
  # Initialize a vector to store the permuted differences
  permuted_delta_fracs <- vector("numeric", n_permutations)
  
  # Perform the permutation test
  library(furrr)
  plan(multisession)
  for (i in 1:n_permutations) {
    # Randomly shuffle the labels
    # when we do not give an argument for size it samples the same size as the input
    # so just shuffling
    shuffled_labels <- sample(test_df[,var])
    # Calculate the difference in the fraction of conserved between the shuffled groups
    delta_frac <- test_df |> 
      mutate(lab = shuffled_labels) %>%
      group_by(lab) |> 
      summarise(N= n(),
                cons_sum = sum(cons),
                frac= cons_sum / n()) |> 
      summarise(delta_frac = diff(frac)) |> 
      pull(delta_frac)
    # Store the permuted difference
    permuted_delta_fracs[i] <- delta_frac
  }
  
  # Calculate the p-value
  p_value <- sum(permuted_delta_fracs >= actual_diff)/n_permutations
  # Print the p-value
  print(p_value)
  
  diffs_df <- data.frame(permuted_diff = permuted_delta_fracs)
  
  p_actual <- actual_df |> 
    ggplot(aes(x = !!sym(var), y = frac)) +
    geom_bar(stat = 'identity', fill = "#000080") +
    ylab('proportion of Smed ATAC peaks in group peak 3') +
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Create the histogram
  histogram <- ggplot(diffs_df , aes(x = permuted_diff)) +
    geom_histogram(bins = 100, fill = "lightblue", color = "black")
  
  p_hist <- histogram +
    geom_vline(aes(xintercept = actual_diff), color = "red", linetype = "dashed", size = 1) +
    theme_minimal() +
    labs(x = "Permuted Differences", y = "Frequency", title = paste0("Permutation Test pvalue: ",
                                                                     format.pval(p_value, digits = 1, eps = 0.001),
         ', ', n_permutations, ' permutations')) +
    annotate("text", x = actual_diff, y = max(ggplot_build(histogram)$data[[1]]$count),
             label = paste("Actual difference =", round(actual_diff, 4)), vjust = -1.5, hjust = 1.5, color = "red") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  pdf(outfile)
 print(
   plot_grid(p_actual, p_hist, nrow = 1,
            rel_widths = c(0.2,1))
 )
  dev.off()
  return(p_value)
}

n_permutations = 10

testing_df <- bi_exp_df |> 
  filter(!is.na(binary_sum)) |> 
  mutate(all_v_promoter = ifelse(element_type == 'putative_promoter', 'YES', 'NO'),
         all_v_enhancer = ifelse(element_type == 'putative_enhancer', 'YES', 'NO'),
         all_v_uncharacterized = ifelse(element_type == 'uncharacterized', 'YES', 'NO'))

permutation_test(testing_df, 'all_v_promoter', n_permutations = n_permutations, 'fig/permutation_test_all_v_promoter.pdf') 



permutation_test(bi_exp_df, 'proximal_CRE', n_permutations = n_permutations, 'fig/permutation_test_proximal_CRE.pdf')

permutation_test(bi_exp_df, 'head', n_permutations = n_permutations, 'fig/permutation_test_head.pdf')
permutation_test(bi_exp_df, 'tail', n_permutations = n_permutations, 'fig/permutation_test_tail.pdf')
permutation_test(bi_exp_df, 'wt',   n_permutations = n_permutations, 'fig/permutation_test_wt.pdf')
permutation_test(bi_exp_df, 'xray', n_permutations = n_permutations, 'fig/permutation_test_xray.pdf')
