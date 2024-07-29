library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)

okabe_ito_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


parse_full_table <- function(tbl, species_name) {
  # TODO: We need to add the chromsize info
  h <- c("busco_id", "status", "sequence", "start", "end", "strand", "score", "length")
  busco <- fread(tbl, fill = TRUE, sep = "\t", skip = 3, header = FALSE)
  names(busco) <- h
  
  # there are interesting modifications of sequence name 
  # This replacement could be avoided by having clean input
  busco2 <- busco %>%
    mutate(sequence = gsub(pattern = "manual_scaffold_", replacement = "chr", sequence))
  busco2$sequence <- gsub(pattern = ":[0-9]+-[0-9]+", replacement = "", busco2$sequence)
  busco2$sp <- species_name
  return(busco2)
}


generate_cosine_curve <- function(x1, y1, x2, y2, n = 100) {
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  x <- seq(x1, x2, length.out = n)
  print('next')
  if(x1 == x2) {
    y = seq(y1, y2, length.out = n)
  } else {
    # Determine the amplitude based on y-values
    amplitude <- (y2 - y1) / 2
    
    # Calculate the necessary phase shift to start at the first y-value
    phase_shift <- pi  # Shift the cosine wave to start from its minimum
    
    # Calculate the y-values along the cosine wave
    y <- y1 + amplitude * cos(((x - x1) / (x2 - x1)) * pi + phase_shift) + amplitude
  }
  return(data.frame(x, y))
}



generate_bezier_curve <- function(x1, y1, x2, y2, n = 100) {
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  x <- seq(x1, x2, length.out = n)
  
  # The distance between x1 and x2
  x_distance <- abs(x2 - x1)
  y_distance <- abs(y2 - y1)
  # use the ratio to place the points
  curvature = min(x_distance / y_distance, 0.3)
  
  
  # Control point height is based on the distance between x1 and x2
  control_height <- (y2 - y1) * 0.3 * sqrt(x_distance)
  
  # Calculate control points
  # Control points are estimated based on the nature of the cosine function
  control1_x <- x1 + (x2 - x1) / 10  # First control point quarter-way
  control2_x <- x1 + 9 * (x2 - x1) / 10  # Second control point three-quarters-way
  
  # this would be in the middle
  # control1_x <- (x1 + x2) / 2
  # control2_x <- (x1 + x2) / 2
  # control1_y <- y1 + control_height
  # control2_y <- y2 + control_height
  
  mid_y <- (y1 + y2) / 2
  amplitude <- (y2 - y1) / 2 # maximum y deviation from mid_y
  
  # we adjust it based on the distance between x1 and x2
  control1_y <- mid_y + amplitude * curvature  # Control point above the line
  control2_y <- mid_y - amplitude * curvature  # Control point below the line
  
  # Generate the Bézier curve using the control points
  t <- seq(0, 1, length.out = n)
  x_bezier <- (1-t)^3 * x1 + 3 * (1-t)^2 * t * control1_x + 3 * (1-t) * t^2 * control2_x + t^3 * x2
  y_bezier <- (1-t)^3 * y1 + 3 * (1-t)^2 * t * control1_y + 3 * (1-t) * t^2 * control2_y + t^3 * y2
  
  curve_data <- data.frame(x = x_bezier, y = y_bezier)
  
  # Return the curve data and control points
  # used for vizualization
  # list(curve_data = curve_data,
  #      control_points = data.frame(x = c(control1_x, control2_x), y = c(control1_y, control2_y)))
  
  return(data.frame(x = x_bezier, y = y_bezier))
}



get_pairwise_cosin <- function(df, sp1, sp2){
  sp1_data <- df %>% filter(plotting_sp == sp1)
  sp2_data <- df %>% filter(plotting_sp == sp2)
  
  # Next, join them by 'busco_id' to get the start and end for each curve
  pairwise_dat <- left_join(sp1_data, sp2_data, by = "busco_id", suffix = c("_start", "_end"))

  # Create an empty data frame to store the curve coordinates
  curve_df <- data.table()
  
  # Loop over each busco_id to generate curves
  for(i in 1:(nrow(pairwise_dat))){
    # Get the start and end points for the current busco_id
    
    # Generate the curve coordinates between the start and end points
    
    curve_points <- generate_bezier_curve(
      pairwise_dat[i,'normalize_start_start'], as.numeric(pairwise_dat[i,'plotting_sp_start']),
      pairwise_dat[i,'normalize_start_end'], as.numeric(pairwise_dat[i,'plotting_sp_end']),
      n = 100)
    
    # curve_points <- generate_cosine_curve(
    #   pairwise_dat[i,'pos_start'], as.numeric(pairwise_dat[i,'plotting_sp_start']),
    #   pairwise_dat[i,'pos_end'], as.numeric(pairwise_dat[i,'plotting_sp_end']),
    #   n = 100)
    
    
    # Add a column for busco_id to the generated points
    curve_points$pair <- paste0(sp1, '_', sp2)
    curve_points$busco_id <- pairwise_dat$busco_id[i]
    curve_points$color <- pairwise_dat$color_start[i]
    # Bind the points to the curve_df
    curve_df <- rbind(curve_df, curve_points)
  }
  
  return(curve_df)
  
  }


################################################################################

## get chromsize

################################################################################

cfiles <- list.files('chromsize/', pattern = '*chromsize')
parse_chromsize <- function(file) {
  df <- read.table(file, header = FALSE)
  names(df) <- c('sequence', 'chromlength')
  df$file <- file
  return(df)
}

lst <- lapply(paste0('chromsize/', cfiles), parse_chromsize) 

all_chromsize <- do.call(rbind, lst) |> 
  mutate(sp = case_when(
    # file == "chromsize/clonorchis_sinensis.PRJNA386618.chromsize" ~ "cloSin",
    # file == "chromsize/hymenolepis_microstoma.PRJEB124.chromsize" ~ "hymMic",
    # file == "chromsize/schistosoma_mansoni.PRJEA36577.chromsize"  ~ "schMan",
    # file == "chromsize/taenia_multiceps.PRJNA307624.chromsize"    ~ "taeMul",
    TRUE ~ gsub('chromsize/(.+).chromsize', '\\1', file)
  ))

res_lst <- list()

genomes = unique(all_chromsize$sp)

res_lst <- lapply(genomes, function(x) {
  parse_full_table(paste0("01_busco/", x, "/busco/run_metazoa_odb10/full_table.tsv"), x)
})

names(res_lst) <- genomes

# add chromsize
res_lst <- lapply(res_lst, left_join, all_chromsize, by = c('sp', 'sequence'))

################################################################################

## Modify the names for plotting

################################################################################

# braLan3 amphioxus
# chromosome names in these assemblies vary quite a bit
# I am here interested in the large contigs that correspond to the chromosomes and I will discard the small contigs
# same goes for haplotigs. For the broad scale analysis we do here they are not needed.

# schPol2
res_lst$schPol2 <- res_lst$schPol2 %>% 
  filter(grepl("chr[1234]$", sequence))
# schMedA2_h1
res_lst$schMedA2_h1 <- res_lst$schMedA2_h1  %>% 
  filter(grepl("Chr", sequence)) %>%
  mutate(sequence = gsub("Chr([0-9])_h.", "chr\\1", sequence))
# schMedA2_h2
res_lst$schMedA2_h2 <- res_lst$schMedA2_h2 %>% 
  filter(grepl("Chr", sequence)) %>%
  mutate(sequence = gsub("Chr([0-9])_h.", "chr\\1", sequence))
# schMedS3_h1
res_lst$schMedS3_h1 <- res_lst$schMedS3_h1 %>% 
  filter(grepl("chr", sequence)) %>%
  mutate(sequence = gsub("chr([0-9])_h.", "chr\\1", sequence))
# schMedS3_h2
res_lst$schMedS3_h2 <- res_lst$schMedS3_h2 %>% 
  filter(grepl("chr", sequence)) %>%
  mutate(sequence = gsub("chr([0-9])_h.", "chr\\1", sequence))

# human GRCh38
res_lst$GRCh38 <- res_lst$GRCh38 %>% 
  filter(grepl("NC_", sequence)) %>%
  mutate(sequence = gsub("NC_[0]+([^.]+)..+", "chr\\1", sequence))
# mouse
res_lst$GRCm39 <- res_lst$GRCm39 %>% 
  filter(grepl("NC_", sequence)) %>%
  mutate(sequence = paste0('chr', as.numeric(gsub("NC_[0]+([^.]+)..+", "\\1", sequence))- 66))

res_lst$GRCm39 <- res_lst$GRCm39

res_lst$GRCm39 
s = 'GRCm39'
all_chromsize |> 
  filter(sp == s) |> 
  arrange(desc(chromlength))

all_chromsize |> 
  group_by(sp) |> 
  summarise(N = n())

# cloSin has 28 chromosomes.trematoda
# I will here take only the numbered ones
res_lst[[7]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst$cloSin <- res_lst$cloSin %>% 
  filter(sequence %in% 1:7) %>%
  mutate(sequence = paste0("chr", sequence))

# hymMic cestoda, six main chromosomes
res_lst[[8]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst$hymMic <- res_lst$hymMic %>% 
  mutate(sequence = gsub(pattern = 'HMN_0(.)_pilon', replacement = 'chr\\1', sequence))

# schMan trematode with six autosomes and a ZW sex chromosome
# Assembly contains haplotypes marked with "H
res_lst[[9]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst$schMan <- res_lst$schMan %>% 
  filter(! (grepl("H0", sequence) | grepl("U0", sequence))) %>%
  mutate(sequence = gsub(pattern = 'SM_V7_', replacement = 'chr', sequence))

# taeMul cestoda, not sure about number of chromosomes but there are only seven that are hit.
res_lst[[10]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst$taeMul <- res_lst$taeMul %>% 
  mutate(sequence = gsub(pattern = 'LG', replacement = 'chr', sequence))


# braLan amphioxous, chromosome scaffolds start with OV and seem to start with OV696686.1
res_lst[[11]] %>%
  group_by(sp, sequence) %>% 
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n()) %>%
  arrange(desc(number_of_hits)) 

res_lst$braLan <- res_lst$braLan %>%
  filter(grepl("OV", sequence)) %>%
  mutate(sequence = gsub(pattern = 'OV(.+).1', replacement = '\\1', sequence)) %>%
  mutate(sequence = as.numeric(sequence) - 696685) %>%
  mutate(sequence = paste0("chr", sequence))

unique(res_lst[[11]]$sequence)

all <- do.call(rbind, res_lst) %>%
  # remove things not called chr from the all to help plotting
  filter(grepl("chr", sequence)) %>%
  # rermove th longer names as well
  filter(! grepl("[0-9]{3}", sequence))

unique(all$sequence)
#grepl("[0-9]{2}", unique(all$sequence))

################################################################################

## output busco coordinates

################################################################################

all %>% 
  select(sp, sequence, start, end, strand, score, length, busco_id, status) %>%
  write.csv(file = "out/busco_locations_2022-busco-syn.csv", row.names = FALSE)


################################################################################

## Plotting Across the Flatworms

################################################################################
sel <- c(
  "schMedS3_h1",
  "schMedS3_h2",
  "schPol2",
  "schNov1",
  "schLug1",
  "cloSin",
  "schMan",
  "hymMic",
  "taeMul",
  "braLan3"
)

parasites <- c("cloSin",
               "schMan",
               "hymMic",
               "taeMul")
step = 0.1
alpha = 0.2
colf <- colorRampPalette(RColorBrewer::brewer.pal(8, name = "Dark2"))
lcol = 'grey50'

################################################################################

## version with schMan as the source of the color

################################################################################
col_source <- all |> 
  filter(sp == 'schMan') |> 
  arrange(sequence, start)

pl <- setNames(okabe_ito_colors, c(paste0('chr', 1:7), 'chrZW'))

col_source$color <- pl[match( col_source$sequence, names(pl))]


all$color <- col_source$color[match(all$busco_id, col_source$busco_id)]

plot_df <- all %>%
  filter(#sequence %in% c("chr1", "chr2", "chr3", "chr4"),
    status == "Complete",
    sp %in% sel) %>%
  mutate(plotting_sp = factor(sp, level = rev(sel)))

ncols <- length(unique(plot_df$busco_id))
col <- colf(ncols)



################################################################################

## normalize the size of the genomes

################################################################################
al = 0.2
fudge_factor = 0.15
fudge_factor = 0.3
new_plot_df <- plot_df |> 
  filter(sp != 'braLan3') |> 
  mutate(
    normalize_start = case_when(
      sequence == "chrZW" ~ 7 + 8 * fudge_factor + round(start/chromlength, 3) ,
      TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 + 
        round(start/chromlength, 3) +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmin = case_when(
      sequence == "chrZW" ~ 7 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmax = case_when(
      sequence == "chrZW" ~ 8 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    )) 

chr_df <- new_plot_df |> 
  select(plotting_sp, normalized_chrmin, normalized_chrmax) |> 
  unique()


new_plot_df <- new_plot_df |> 
  filter(!is.na(color)) |> 
  filter(sp != 'braLan3') 

sc_complete <- new_plot_df |> 
  group_by(busco_id) |> 
  summarise(N = n()) |> 
  filter(N == 9) |> 
  pull(busco_id)

# UNIQE filter
new_plot_df <- new_plot_df |>
  filter(busco_id %in% sc_complete)

pdf('fig/busco_syn_flatworms_normalized_schMan_colors_grey_v2.pdf', width = 8, height = 10)

new_plot_df %>% 

  ggplot(aes(x=normalize_start, y = plotting_sp, color = color)) +
  geom_point(size = 0.5) +
  #geom_vline(xintercept = 1:8) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h1", "schMedS3_h2")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h2", "schPol2")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schPol2", "schNov1")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schNov1", "schLug1")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schLug1", "cloSin")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("cloSin", "schMan")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMan", "hymMic")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("hymMic", "taeMul")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_line(data = filter(new_plot_df,  sp %in% c("taeMul", "braLan3")),
            aes(group = busco_id) ,alpha = al, color =  lcol) +
  geom_segment(data = chr_df,
               aes(y = as.numeric(plotting_sp) -1, # some issue with the factor levels...
                   yend= as.numeric(plotting_sp) -1,
                   x = normalized_chrmin,
                   xend=normalized_chrmax),
               color = 'grey30', size = 4, lineend='round') +
  geom_segment(aes(y = as.numeric(plotting_sp) - 0.2 - 1,
                   yend= as.numeric(plotting_sp) + 0.2 - 1,
                   xend=normalize_start)) +
  scale_color_manual(values = c( "#E69F00","#F0E442", "#D55E00", "#0072B2", "#CC79A7","#000000","#009E73",  "#56B4E9")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) 

dev.off()


#################################################################################

## with curves
links <- list(
  get_pairwise_cosin(new_plot_df, "schMedS3_h1", "schMedS3_h2"),
  get_pairwise_cosin(new_plot_df, "schMedS3_h2", "schPol2"),
  get_pairwise_cosin(new_plot_df, "schPol2", "schNov1"),
  get_pairwise_cosin(new_plot_df, "schNov1", "schLug1"),
  
  get_pairwise_cosin(new_plot_df, "schLug1", "cloSin"),
  get_pairwise_cosin(new_plot_df, "cloSin", "schMan"),
  get_pairwise_cosin(new_plot_df, "schMan", "hymMic"),
  get_pairwise_cosin(new_plot_df, "hymMic", "taeMul")
)

al = 0.8
sz = 0.5
png('fig/busco_syn_flatworms_bezier_schMan_colors.png', width = 10, height = 10,
    units = 'in', res = 900)

new_plot_df %>% 
  ggplot(aes(x=normalize_start, y = plotting_sp, color = color)) +
  geom_point(size = 0.5) +
  #geom_vline(xintercept = 1:8) +
  geom_line(data = links[[1]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_line(data = links[[2]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_line(data = links[[3]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_line(data = links[[4]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_line(data = links[[5]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_line(data = links[[6]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_line(data = links[[7]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_line(data = links[[8]], aes(x, y - 1 , group = busco_id), alpha = al, size = sz) +
  geom_segment(data = chr_df,
               aes(y = as.numeric(plotting_sp) -1, # some issue with the factor levels...
                   yend= as.numeric(plotting_sp) -1,
                   x = normalized_chrmin,
                   xend=normalized_chrmax),
               color = 'grey30', size = 5, lineend='round') +
  geom_segment(aes(y = as.numeric(plotting_sp) - 0.07 - 1,
                   yend= as.numeric(plotting_sp) + 0.07 - 1,
                   xend=normalize_start)) +
  scale_color_manual(values = c( "#E69F00","#F0E442", "#D55E00", "#0072B2", "#CC79A7","#000000","#009E73",  "#56B4E9")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) 


dev.off()

######

################################################################################


pdf('fig/busco_syn_flatworms_normalized_schMan_colors.pdf', width = 8, height = 10)

new_plot_df %>% 
  ggplot(aes(x=normalize_start, y = plotting_sp, color = color)) +
  geom_point(size = 0.5) +
  #geom_vline(xintercept = 1:8) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h1", "schMedS3_h2")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h2", "schPol2")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schPol2", "schNov1")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schNov1", "schLug1")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schLug1", "cloSin")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("cloSin", "schMan")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMan", "hymMic")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("hymMic", "taeMul")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("taeMul", "braLan3")),
            aes(group = busco_id) ,alpha = al) +
  geom_segment(data = chr_df,
               aes(y = as.numeric(plotting_sp) -1, # some issue with the factor levels...
                   yend= as.numeric(plotting_sp) -1,
                   x = normalized_chrmin,
                   xend=normalized_chrmax),
               color = 'grey30', size = 4, lineend='round') +
  geom_segment(aes(y = as.numeric(plotting_sp) - 0.05 - 1,
                   yend= as.numeric(plotting_sp) + 0.05 - 1,
                   xend=normalize_start)) +
  scale_color_manual(values = c( "#E69F00","#F0E442", "#D55E00", "#0072B2", "#CC79A7","#000000","#009E73",  "#56B4E9")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) 

dev.off()








################################################################################

## Below is exploratory plotting

################################################################################

################################################################################

## Basic summary stats

################################################################################

all %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n()) %>%
  View()

################################################################################

## Against Human

################################################################################
a <- do.call(rbind, res_lst)


################################################################################
## color setup
col_source <- a |> 
  filter(grepl('chr', sequence), sp == "braLan3") |> 
  arrange(sequence, start)

n_chromosomes <- length(unique(col_source$sequence))
color_palette <- hcl.colors(n_chromosomes, "viridis")
chromosome_names <- paste0('chr', 1:n_chromosomes)
pl <- setNames(color_palette, chromosome_names)

col_source$color <- pl[match( col_source$sequence, names(pl))]


a$color <- col_source$color[match(a$busco_id, col_source$busco_id)]


sel <- c(
  "GRCh38"

)

sel <- c(
  "braLan3",
  "GRCh38",
  "GRCm39"
)

plot_df <- a %>%
  filter(#sequence %in% c("chr1", "chr2", "chr3", "chr4"),
    status == "Complete",
    sp %in% sel) %>%
  mutate(plotting_sp = factor(sp, level = rev(sel)))


fudge_factor = 0.6
new_plot_df <- plot_df |> 
  #filter(sp != 'braLan3') |> 
  mutate(
    normalize_start = case_when(
      sequence == "chrZW" ~ 7 + 8 * fudge_factor + round(start/chromlength, 3) ,
      TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 + 
        round(start/chromlength, 3) +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmin = case_when(
      sequence == "chrZW" ~ 7 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmax = case_when(
      sequence == "chrZW" ~ 8 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    )) 


chr_df <- new_plot_df |> 
  select(plotting_sp, normalized_chrmin, normalized_chrmax) |> 
  unique()

new_plot_df %>% 
  ggplot(aes(x=normalize_start, y = plotting_sp, color = color)) +
  geom_point(size = 0.5) +
  #geom_vline(xintercept = 1:8) +
  geom_line(data = filter(new_plot_df,  sp %in% c("GRCh38",
                                                  "braLan3")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("GRCh38",
                                                  "GRCm39")),
            aes(group = busco_id) ,alpha = al) +
  geom_segment(data = chr_df,
               aes(y = as.numeric(plotting_sp), # some issue with the factor levels...
                   yend= as.numeric(plotting_sp),
                   x = normalized_chrmin,
                   xend=normalized_chrmax),
               color = 'grey30', size = 4, lineend='round') +
  geom_segment(aes(y = as.numeric(plotting_sp) - 0.05,
                   yend= as.numeric(plotting_sp) + 0.05,
                   xend=normalize_start)) +
  #scale_color_brewer(palette = 'Dark2') +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) +
  scale_y_discrete(labels = rev(c('Amophioxus', 'Human', 'Mouse')))




sel <- c(
  "braLan3",
  "GRCh38"
)

plot_df <- a %>%
  filter(#sequence %in% c("chr1", "chr2", "chr3", "chr4"),
    # status == "Complete",
    sp %in% sel) %>%
  mutate(plotting_sp = factor(sp, level = rev(sel))) |> 
  group_by(busco_id) |> 
  mutate(running_number = row_number())

new_plot_df <- plot_df |> 
  #filter(sp != 'braLan3') |> 
  mutate(
    normalize_start = case_when(
      sequence == "chrZW" ~ 7 + 8 * fudge_factor + round(start/chromlength, 3) ,
      TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 + 
        round(start/chromlength, 3) +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmin = case_when(
      sequence == "chrZW" ~ 7 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmax = case_when(
      sequence == "chrZW" ~ 8 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    )) 

g1 <- plot_df |> 
  filter(sp == 'braLan3',
         grepl('chr', sequence)) |> 
  select(busco_id, sequence, start, end, color, chromlength) |> 
  rename(g1_sequence = sequence,
         g1_start = start,
         g1_end = end,
         g1_chromlength = chromlength)

chromosomes <- (unique(g1$g1_sequence))
sorted_chromosomes <- chromosomes[order(as.numeric(gsub("chr", "", chromosomes)))]
g1$g1_sequence <- factor(g1$g1_sequence, levels = sorted_chromosomes)

g2 <- plot_df |> 
  filter(sp == "GRCh38",
         grepl('chr', sequence)) |> 
  select(busco_id, sequence, start, end, color, chromlength) |> 
  rename(g2_sequence = sequence,
         g2_start = start,
         g2_end = end,
         g2_chromlength = chromlength)

chromosomes <- (unique(g2$g2_sequence))
sorted_chromosomes <- chromosomes[order(as.numeric(gsub("chr", "", chromosomes)))]
g2$g2_sequence <- factor(g2$g2_sequence, levels = sorted_chromosomes)


g <- full_join(g1, g2, by = c('busco_id',  'color'), multiple = 'all')

g |> 
  na.omit() |> 
  ggplot(aes(x = g2_start, y = g1_start, color = color)) +
  geom_point() + 
  facet_grid(rows = vars(g2_sequence), cols  = vars(g1_sequence),
             scales = 'free', space = 'free', margins = FALSE) +
  theme_bw() +
  theme(text = element_text(size = 8)) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    # axis.text.y = element_text(size = 10, angle = 15),
    # axis.text.x = element_text(size = 10, angle = 15),
    legend.position = 'none',
    strip.background = element_rect(fill = NA),
    text = element_text(size = 12),
    panel.spacing = unit(0, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab('braLan3') +
  ylab('GRCh38')



plot_df <- new_plot_df |> 
  filter(grepl('chr', sequence)) |> 
  select(sequence, start, sp, busco_id, color) |> 
  pivot_wider(id_cols = c(busco_id, color), names_from =  sp,
              values_from = c(start, sequence)) |> 
  na.omit()


new_plot_df %>%
  dplyr::group_by(busco_id, color, sp) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 









plot_df  %>% 
  filter(sp != 'braLan3') |> 
  ggplot(aes(x=normalize_start, y = plotting_sp, color = color)) +
  geom_point(size = 0.5) +
  #geom_vline(xintercept = 1:8) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h1", "schMedS3_h2")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h2", "schPol2")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schPol2", "schNov1")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schNov1", "schLug1")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schLug1", "cloSin")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("cloSin", "schMan")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMan", "hymMic")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("hymMic", "taeMul")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("taeMul", "braLan3")),
            aes(group = busco_id) ,alpha = al) +
  geom_segment(data = chr_df,
               aes(y = as.numeric(plotting_sp) -1, # some issue with the factor levels...
                   yend= as.numeric(plotting_sp) -1,
                   x = normalized_chrmin,
                   xend=normalized_chrmax),
               color = 'grey30', size = 4, lineend='round') +
  geom_segment(aes(y = as.numeric(plotting_sp) - 0.05 - 1,
                   yend= as.numeric(plotting_sp) + 0.05 - 1,
                   xend=normalize_start)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) 

################################################################################

## only S3_h1

################################################################################



################################################################################

## dotplot intermezzo

################################################################################



unique(df$sequence)
dd <- df |> 
  filter(sp %in% c('GRCh38')) 

unique(dd$sequence)


plot_df <- all %>%
  filter(#sequence %in% c("chr1", "chr2", "chr3", "chr4"),
    status == "Complete",
    sp %in% c('GRCh38', 'GRCm39')) %>%
  mutate(plotting_sp = factor(sp, level = rev(sel)))




################################################################################

## color the BUSCO genes by chromosome of S3h1

################################################################################

col_source <- all |> 
  filter(sp == 'schMedS3_h1') |> 
  arrange(sequence, start)

pl <- setNames(RColorBrewer::brewer.pal(n=4, 'Dark2'), paste0('chr', 1:4))

col_source$color <- pl[match( col_source$sequence, names(pl))]


all$color <- col_source$color[match(all$busco_id, col_source$busco_id)]

plot_df <- all %>%
  filter(#sequence %in% c("chr1", "chr2", "chr3", "chr4"),
    status == "Complete",
    sp %in% sel) %>%
  mutate(plotting_sp = factor(sp, level = rev(sel)))

ncols <- length(unique(plot_df$busco_id))
col <- colf(ncols)





unique(plot_df$sp)

################################################################################

## normalize the size of the genomes

################################################################################
al = 0.2
fudge_factor = 0.15
new_plot_df <- plot_df |> 
  filter(sp != 'braLan3') |> 
  mutate(
    normalize_start = case_when(
    sequence == "chrZW" ~ 7 + 8 * fudge_factor + round(start/chromlength, 3) ,
    TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 + 
      round(start/chromlength, 3) +
      (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmin = case_when(
      sequence == "chrZW" ~ 7 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) - 1 +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    ),
    normalized_chrmax = case_when(
      sequence == "chrZW" ~ 8 + 8 * fudge_factor,
      TRUE ~ as.numeric(gsub('chr','',sequence)) +
        (as.numeric(gsub('chr','',sequence)) * fudge_factor)
    )) 

chr_df <- new_plot_df |> 
  select(plotting_sp, normalized_chrmin, normalized_chrmax) |> 
  unique()

pdf('fig/busco_syn_flatworms_normalized.pdf', width = 8, height = 10)
new_plot_df %>% 
  filter(sp != 'braLan3') |> 
  ggplot(aes(x=normalize_start, y = plotting_sp, color = color)) +
  geom_point(size = 0.5) +
  #geom_vline(xintercept = 1:8) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h1", "schMedS3_h2")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h2", "schPol2")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schPol2", "schNov1")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schNov1", "schLug1")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schLug1", "cloSin")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("cloSin", "schMan")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("schMan", "hymMic")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("hymMic", "taeMul")),
            aes(group = busco_id) ,alpha = al) +
  geom_line(data = filter(new_plot_df,  sp %in% c("taeMul", "braLan3")),
            aes(group = busco_id) ,alpha = al) +
  geom_segment(data = chr_df,
               aes(y = as.numeric(plotting_sp) -1, # some issue with the factor levels...
                   yend= as.numeric(plotting_sp) -1,
                   x = normalized_chrmin,
                   xend=normalized_chrmax),
  color = 'grey30', size = 4, lineend='round') +
  geom_segment(aes(y = as.numeric(plotting_sp) - 0.05 - 1,
                   yend= as.numeric(plotting_sp) + 0.05 - 1,
                   xend=normalize_start)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) 

 dev.off()
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
 lcol = 'grey60'
 
 pdf('fig/busco_syn_flatworms_normalized_grey.pdf', width = 8, height = 10)
 new_plot_df %>% 
   filter(sp != 'braLan3') |> 
   ggplot(aes(x=normalize_start, y = plotting_sp, color = color)) +
   geom_point(size = 0.5) +
   #geom_vline(xintercept = 1:8) +
   geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h1", "schMedS3_h2")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("schMedS3_h2", "schPol2")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("schPol2", "schNov1")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("schNov1", "schLug1")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("schLug1", "cloSin")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("cloSin", "schMan")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("schMan", "hymMic")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("hymMic", "taeMul")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_line(data = filter(new_plot_df,  sp %in% c("taeMul", "braLan3")),
             aes(group = busco_id) ,alpha = al, color =  lcol) +
   geom_segment(data = chr_df,
                aes(y = as.numeric(plotting_sp) -1, # some issue with the factor levels...
                    yend= as.numeric(plotting_sp) -1,
                    x = normalized_chrmin,
                    xend=normalized_chrmax),
                color = 'grey30', size = 4, lineend='round') +
   geom_segment(aes(y = as.numeric(plotting_sp) - 0.05 - 1,
                    yend= as.numeric(plotting_sp) + 0.05 - 1,
                    xend=normalize_start)) +
   scale_color_brewer(palette = 'Dark2') +
   theme_minimal() +
   theme(legend.position = "none",
         axis.title = element_blank(),
         axis.text.x = element_blank(),
         line = element_blank()) 
 
 dev.off() 
 





################################################################################
## TESTING
# filter to only contain those that are found in all


rip_df <- new_plot_df %>% 
  filter(sp != 'braLan3') |> 
  #filter(sp %in% c("schMedS3_h1", "schNov1","schLug1", "schMan")) |> 
  mutate(plotting_sp = factor(plotting_sp)) |> 
  mutate(num_sp = as.numeric(plotting_sp)) 

complete_buscos <- rip_df |> 
  group_by(busco_id) |> 
  summarise(N = n()) |> 
  filter(N == max(N)) |> 
  pull(busco_id)

f <- c('131349at33208')
spline_s = 50
spline_m = "natural"
# https://stackoverflow.com/questions/54898700/connecting-points-with-curved-line-on-ggplot-for-a-categorical-variable-on-the
xvar = 'plotting_sp'
yvar = 'normalize_start'

rip_df |> 
  filter(busco_id %in% complete_buscos) |> 
  ggplot(aes(y=!!sym(yvar), x = !!sym(xvar), color = color)) +
  #geom_point(size = 0.5) +
  geom_line(
    data = . %>%
      mutate(plotting_sp = as.numeric(!!sym(xvar))) %>%
      group_by(busco_id, color) %>%
      # increase n for smoother line; can also try the other methods listed under
      # ?spline, though I find "natural" looks better than some of the rest: "fmm"'s 
      # curves are rather drastic, & "periodic" doesn't touch all the points
      summarise(x1 = list(spline(!!sym(xvar), !!sym(yvar), n = spline_s, method = spline_m)[["x"]]),
                y1 = list(spline(!!sym(xvar), !!sym(yvar), n = spline_s, method = spline_m)[["y"]])) %>%
      tidyr::unnest(cols = c(x1, y1)),
    aes(x = x1, y = y1, color = color, group = busco_id),
    alpha = al) +
  # geom_segment(data = chr_df,
  #              aes(x = as.numeric(plotting_sp) -1, # some issue with the factor levels...
  #                  xend= as.numeric(plotting_sp) -1,
  #                  y = normalized_chrmin,
  #                  yend=normalized_chrmax),
  #              color = 'grey30', size = 4, lineend='round') +
  geom_segment(aes(x = as.numeric(plotting_sp) - 0.05 ,
                   xend= as.numeric(plotting_sp) + 0.05 ,
                   yend=normalize_start)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) +
  coord_flip()

# bit done testing