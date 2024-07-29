library(tidyverse)
library(writexl)

format_posthoc <- function(path, out) {
  df <- read.csv(path)
  df |> 
    mutate(Value = gsub('p value', 'pval', Value)) |> 
    pivot_wider(id_cols =  Dimension, names_from = Value, values_from = -c(X, Dimension, Value)) |>
    # Round Residuals
    mutate(across(ends_with('Residuals'), function(x){round(as.numeric(x), 1)})) %>%
    # Remove asterisks from pvals
    mutate(across(ends_with('pvals'), ~ gsub("\\*", "", .))) %>%
    # Convert pvals columns to numeric
    mutate(across(ends_with('pvals'), function(x){
      format.pval(as.numeric(x), eps = 1e-16, digits = 3)
    })) |> 
    write_xlsx(out)
}

format_posthoc('./peak_annotation/Schmidtea_annotation_by_species_chi_posthoc.csv',
               './peak_annotation/Schmidtea_annotation_by_species_chi_posthoc.xlsx')

format_posthoc('./peak_annotation/schMedS3h1_distribution_by_element_type_chi_posthoc.csv',
               './peak_annotation/schMedS3h1_distribution_by_element_type_chi_posthoc.xlsx')

# conservation 
format_posthoc('./conserved_elements/chi_posthoc.csv',
               './conserved_elements/chi_posthoc.xlsx')
