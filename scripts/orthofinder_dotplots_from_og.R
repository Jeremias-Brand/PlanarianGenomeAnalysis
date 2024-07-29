library(tidyverse)
library(ggplot2)
library(ape)
library(RColorBrewer)
library(data.table)

bed <- fread('orthofinder/combBed.txt')
mb_contigs <- bed |> 
  group_by(genome, chr) |> 
  summarise(maxLen = max(c(start, end))) |> 
  filter(maxLen > 10e6)


make_dotplot <- function(bed, target, query, scales = 'free', space = 'free', palette = 'black') {
  mb_contigs <- bed |> 
    group_by(genome, chr) |> 
    summarise(maxLen = max(c(start, end))) |> 
    filter(maxLen > 10e6)
  
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
    summarise(N = n()) |> 
    filter(N <= 1) |> 
    pull(globOG)
  
  

  joined_sc_df <- joined_df |> 
    filter(globOG %in% sc_ogs) |> 
    na.omit() |> 
    mutate(chr_query = factor(chr_query))
  
  joined_sc_df$chr_query <- factor(joined_sc_df$chr_query, levels = rev(levels(joined_sc_df$chr_query)))
  
  # set3_palette <- brewer.pal(12, "Set3")
  # additional_colors <- c("#3D3D3D", "#B30000", "#006E00", "#0000B3", "#B35900", "#800080", "#990099")
  # custom_palette <- c(set3_palette, additional_colors)
  if (length(palette) == 1) {
    p <- ggplot(joined_sc_df, aes(x = start_target, y = start_query)) +
      geom_point(size = 0.2, alpha = 1, color = palette)
  } else {
    p <- ggplot(joined_sc_df, aes(x = start_target, y = start_query, color = factor(chr_target))) +
      geom_point(size = 0.2, alpha = 1)
  }
  
   # p <- p +
   #  facet_grid(chr_query ~ chr_target, scales = scales, space = space) +
   #  theme_minimal() +
   #  scale_x_continuous(expand = c(0,0)) +
   #  scale_y_continuous(expand = c(0,0)) +
   #  labs(title = paste0('Target (x): ', target, ' Query (y): ', query),
   #       x = target,
   #       y = query) +
   #  scale_color_manual(values = custom_palette, name = "Chromosome (Query)") +
   #  theme(legend.position = "none",
   #        axis.text = element_blank(),
   #        panel.grid = element_blank(),
   #        axis.ticks = element_blank(),
   #        panel.spacing.x = unit(0.0, "lines"),  # Reduce the distance between facets horizontally
   #        panel.spacing.y = unit(0.0, "lines"),
   #        panel.border = element_rect(color = "grey10", fill = NA), 
   #        strip.background = element_rect(fill = "white"),  # Set facet label background to white
   #        strip.text = element_text(size = 8) )
   
   p <- p +
     facet_grid(chr_query ~ chr_target, scales = scales, space = space) +
     theme_minimal() +
     scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), breaks = scales::breaks_pretty(n = 5), expand = c(0,0)) +
     scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), breaks = scales::breaks_pretty(n = 5), expand = c(0,0)) +
     labs(x = target,
          y = query) +
     scale_color_manual(values = custom_palette, name = "Chromosome (Query)") +
     theme(legend.position = "none",
           axis.text = element_blank(),
           panel.grid = element_blank(),
           axis.ticks = element_blank(),
           panel.spacing.x = unit(0.0, "lines"),  # Reduce the distance between facets horizontally
           panel.spacing.y = unit(0.0, "lines"),
           panel.border = element_rect(color = "grey10", fill = NA), 
           strip.background = element_rect(fill = "white"),  # Set facet label background to white
           strip.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
           strip.text.y = element_text(size = 8, angle = 0, vjust = 0.5)
     )

   
  
  return(p)
  
}


 make_dotplot_main <- function(bed, target, query, scales = 'free', space = 'free', palette = 'black') {
  mb_contigs <- bed |>
    group_by(genome, chr) |>
    summarise(maxLen = max(c(start, end))) |>
    filter(maxLen > 10e6)
  
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
    summarise(N = n()) |> 
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
  

  if (length(palette) == 1) {
    p <- ggplot(joined_sc_df, aes(x = start_target, y = start_query)) +
      geom_point(size = 0.3, alpha = 1, color = palette)
  } else {
    p <- ggplot(joined_sc_df, aes(x = start_target, y = start_query, color = factor(chr_target))) +
      geom_point(size = 0.3, alpha = 1)
  }
  
  p <- p +
    facet_grid(chr_query ~ chr_target, scales = scales, space = space) +
    theme_minimal() +
    scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), breaks = scales::breaks_pretty(n = 5), expand = c(0,0)) +
    scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), breaks = scales::breaks_pretty(n = 5), expand = c(0,0)) +
    # scale_x_continuous(expand = c(0,0)) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(x = '',
         y = query) +
    scale_color_manual(values = custom_palette, name = "Chromosome (Query)") +
    theme(legend.position = "none",
          axis.text = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing.x = unit(0.0, "lines"),  # Reduce the distance between facets horizontally
          panel.spacing.y = unit(0.0, "lines"),
          panel.border = element_rect(color = "grey10", fill = NA), 
          strip.background = element_rect(fill = "white"),  # Set facet label background to white
          strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 8, angle = 0, vjust = 0.5)
          )
  
  return(p)
 }
 

custom_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
                    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#BEBADA", "black", brewer.pal(12, "Set3"))



gen <- unique(bed$genome)
targets <- c("BraLan", "schMedS3h1", "schMan")

plot_list <- list()

for (t in targets) {
  for (q in gen) {
    plot_list[[paste(t, q, sep = "_")]] <- make_dotplot(bed, t, q)
  }
}

png_out = 'orthofinder/og_dotplots/'


for (plot_name in names(plot_list)) {
  plot <- plot_list[[plot_name]]
  png(paste0(png_out,plot_name,'.png'),res = 300, units = 'in', width = 6, height = 6)
  print(plot)
  dev.off()
}

# Create a PDF file
pdf("orthofinder/orthofinder_dotplots_grayscale.pdf", width = 10, height = 8)

for (plot_name in names(plot_list)) {
  plot <- plot_list[[plot_name]]
  print(plot)
}

dev.off()





