library(tidyverse)
# get the file names that match the pattern
fs <- list.files('02_merqury/', pattern = '*completeness.stats')

# loop through the files
for (f in fs) {
  df <- read.table(paste0('02_merqury/', f), header = FALSE)
  name <- sub("_.*", "", f)
  df$name <- name
  df_list[[length(df_list) + 1]] <- df
}

result_df <- do.call(rbind, df_list)
names(result_df) <- c('assembly', 'group', 'all_kmer', 'solid_kmer_in_assembly', 'completeness', 'name')


qv_list <- list()

# get the file names that match the pattern
# I only want the overall results not split by chromosome
fs <- list.files('02_merqury/', pattern = 'SRR[^.]+.qv')

for (f in fs) {
  df <- read.table(paste0('02_merqury/', f), header = FALSE)
  name <- sub("_.*", "", f)
  df$name <- name
  qv_list[[length(qv_list) + 1]] <- df
}

# combine all the dataframes in the list into a single dataframe
qv_df <- do.call(rbind, qv_list)
names(qv_df) <- c('assembly', 'kmer_assembly_only', 'kmerassembly_reads', 'qv', 'error_rate', 'name')

df <- full_join(result_df, qv_df)


qv_p <- df |> 
  filter(!grepl('A2_h', assembly) & !grepl('GCA', assembly)) |> 
ggplot(data = , 
       aes(x = completeness, y = qv, color = assembly, shape = name)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(19, 17, 15, 21)) + # customize shape according to your need
  labs(x = "Completeness", y = "Quality Value", color = "Assembly", shape = "Illumina data") +
  theme_bw() +
  scale_color_brewer(palette = 'Dark2')

pdf('fig/QV_schmidtea_sexual_assemblies.pdf', width = 10)

col = c("#E7298A", "#B5216D", "#A7A1FF", "#7570B3" , "#6054FF")

df |> 
  filter(!grepl('A2_h', assembly) & !grepl('GCA', assembly)) |> 
  ggplot(data = , 
         aes(x = completeness, y = qv, color = assembly, shape = name)) +
  geom_point(size = 8) +
  scale_shape_manual(values = c(19, 17, 15, 21)) + # customize shape according to your need
  labs(x = "Completeness", y = "Quality Value", color = "Assembly", shape = "Illumina data") +
  theme_bw() +
  scale_color_manual(values = col) +
  theme(text = element_text(size=21))
dev.off()
