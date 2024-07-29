library(tidyverse)
library(ggplot2)
library(ape)
library(RColorBrewer)
library(data.table)
library(furrr)
library(parallel)

################################################################################

## functions

################################################################################

filter_bed <- function(bed, target, query) {
  pair = c(target, query)
  
  filtered_bed <- bed |> 
    select(chr, start, end, globOG, genome) |> 
    filter(genome %in% pair, chr %in% mb_contigs$chr) |> 
    group_by(genome) 
  
  unique_genomes <- filtered_bed %>% distinct(genome) %>% pull(genome)
  
  filtered_bed_list <- filtered_bed |> 
    group_split() 
  
  names(filtered_bed_list) <- unique_genomes
  
  # The first dataframe is for the 'target' genome, and the second is for the 'query' genome
  target_df <- filtered_bed_list[[match(target, names(filtered_bed_list))]]
  if (target == query) {
    query_df = target_df
  } else {
    query_df <- filtered_bed_list[[match(query, names(filtered_bed_list))]]
  }
  
  # Perform inner join on the 'globOG' variable
  joined_df <- target_df %>%
    left_join(query_df, by = "globOG", suffix = c("_target", "_query"), multiple = 'all') |> 
    filter(!grepl('NoOG', globOG))
  
  
  sc_ogs <- joined_df |> 
    group_by(globOG) |> 
    summarise(N = dplyr::n()) |> 
    filter(N <= 1) |> 
    pull(globOG)
  
  joined_sc_df <- joined_df |> 
    filter(globOG %in% sc_ogs) |> 
    na.omit() |> 
    mutate(chr_query = factor(chr_query))
  
  joined_sc_df$chr_query <- factor(joined_sc_df$chr_query, levels = rev(levels(joined_sc_df$chr_query)))
  
  ################################################################################
  ## ordering
  # Extract numeric part from chr_query
  joined_sc_df$chr_query_num <- as.numeric(gsub("chr", "", joined_sc_df$chr_query))
  # Reorder chr_query based on the numeric order
  joined_sc_df$chr_query <- factor(joined_sc_df$chr_query, levels = rev(unique(joined_sc_df$chr_query[order(joined_sc_df$chr_query_num)])))
  # Extract numeric part from chr_target
  joined_sc_df$chr_target_num <- as.numeric(gsub("chr", "", joined_sc_df$chr_target))
  # Reorder chr_target based on the numeric order
  joined_sc_df$chr_target <- factor(joined_sc_df$chr_target, levels = unique(joined_sc_df$chr_target[order(joined_sc_df$chr_target_num)]))
  
  return(joined_sc_df)
}

permutation_chi2 <- function(joined_sc_df, workers = 4, permutations = 10000) {
  # Compute the observed contingency table

  observed_table <- table(joined_sc_df$chr_target, joined_sc_df$chr_query)
  
  # Calculate row and column totals
  row_totals <- margin.table(observed_table, 1)
  col_totals <- margin.table(observed_table, 2)
  
  # Calculate the grand total
  grand_total <- sum(observed_table)
  
  # Compute the expected contingency table
  expected <- outer(row_totals, col_totals, "*") / grand_total
  
  # Rest of the permutation test as provided before:
  
  # Define the number of permutations
  n_permutations <- permutations
  
  # Define the observed test statistic (chi-square statistic)
  observed_statistic <- sum((observed_table - expected)^2 / expected)
  
  # Initialize a vector to store the permuted test statistics
  permuted_statistics <- numeric(n_permutations)
  
  
  
  plan(multisession, workers = workers)
  
  # Initialize a vector to store the permuted test statistics
  print(n_permutations)
  permuted_statistics <- future_map_dbl(1:n_permutations, function(i) {
    # Permute the chr_query column
    permuted_chr_query <- sample(joined_sc_df$chr_query)
    
    # Compute the permuted contingency table
    permuted_table <- table(joined_sc_df$chr_target, permuted_chr_query)
    
    # Compute the permuted test statistic (chi-square statistic)
    permuted_statistic <- sum((permuted_table - expected)^2 / expected)
    
    # Return the permuted statistic
    return(permuted_statistic)
  }, .progress = TRUE)
  
  # Compute the p-value as the proportion of permuted test statistics that are greater than the observed test statistic
  p_value <- mean(permuted_statistics >= observed_statistic)
  if (p_value == 0) {
    p_value = paste0('< ', 1/n_permutations)
  }
  
  return(p_value)
}

# Define a function to perform the chi-square test and calculate Cramer's V
perform_chi_square <- function(df, query, target) {
  # Create a contingency table
  chromosome_table <- table(df$chr_target, df$chr_query)
  
  # Perform the chi-square test
  chi2_test_result <- chisq.test(chromosome_table)
  
  # Compute Cramer's V
  n <- sum(chromosome_table)  # total number of observations
  k <- min(dim(chromosome_table))  # number of rows or columns (whichever is smaller)
  V <- sqrt(chi2_test_result$statistic / (n*(k-1)))
  
  # Return a list with all the results
  data.frame(
    target = target,
    query = query, 
    N = n,
    df = chi2_test_result$parameter,
    chisq = round(chi2_test_result$statistic, 1),
    V = round(unname(V),2),
    pval = format.pval(chi2_test_result$p.value, digits = 3, eps = 1e-16)
  )
}


bed <- fread('orthofinder/combBed.txt')
mb_contigs <- bed |> 
  group_by(genome, chr) |> 
  summarise(maxLen = max(c(start, end))) |> 
  filter(maxLen > 10e6)

targets <- c("BraLan", "NemVec", "schMedS3h1", "schMan")

queries <- c("BraLan", "NemVec", "cloSin","echMul","hymMic","schMan" ,"taeMul", 
             "schMedS3h1", "schMedS3h2", "schPol2", "schNov1", "schLug1")

################################################################################

## statistical testing

################################################################################

################################################################################
## vanilla chi-square test

res <- list()
for (target in targets) {
  for (query in queries) {
    #. stats
    f_bed <- filter_bed(bed, target = target, query = query)
    res[[paste0(target, '_', query)]] <- perform_chi_square(f_bed, target = target, query = query)
    
    # Convert the table to a data frame in a long format
    chromosome_table <- table(f_bed$chr_target, f_bed$chr_query)
    chromosome_df <- as.data.frame(chromosome_table)
    names(chromosome_df) <- c("chr_target", "chr_query", "Frequency")
    
    # Plot the heatmap
    pdf(paste0('orthofinder/synteny_chisq_', target, query, '.pdf'))
    print(
      ggplot(chromosome_df, aes(chr_target, chr_query, fill = Frequency)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(x = "Chr Target", y = "Chr Query", fill = "Frequency", 
           title = "Heatmap of Ortholog Distribution Across Chromosomes")
    )
    dev.off()

    
  }
}

chisq_res <- do.call(rbind, res)

################################################################################
## Simulation based p-values
permutations = 100000
pvals <- list()

for (target in targets) {
  for (query in queries) {
    print(target)
    print(query)
    f_bed <- filter_bed(bed, target = target, query = query)
    pvals[[paste0(target, '_', query)]] <- permutation_chi2(
      joined_sc_df = f_bed, workers = 7, permutations = permutations
      ) 
  }
}

perm_res <- do.call(rbind, pvals) |> 
  as.data.frame() |> 
  rownames_to_column() |> 
  rowwise() |> 
  mutate(target = strsplit(rowname, split = '_')[[1]][[1]],
         query = strsplit(rowname, split = '_')[[1]][[2]]) |> 
  select(-rowname) |> 
  dplyr::rename(perm_pval = V1)

perm_res$permutations = permutations
chisq_res <- left_join(chisq_res, perm_res)

write.table(chisq_res, file = 'orthofinder/synteny_chisq.tsv', 
            quote = FALSE, row.names = FALSE, sep = '\t')











