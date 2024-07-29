library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)


parse_full_table <- function(tbl, species_name) {
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
    file == "chromsize/clonorchis_sinensis.PRJNA386618.fa.chromsize" ~ "cloSin",
    file == "chromsize/hymenolepis_microstoma.PRJEB124.fa.chromsize" ~ "hymMic",
    file == "chromsize/schistosoma_mansoni.PRJEA36577.fa.chromsize"  ~ "schMan",
    file == "chromsize/taenia_multiceps.PRJNA307624.fa.chromsize"    ~ "taeMul",
    file == "chromsize/schMedA2_h1.fa.chromsize" ~ "schMedA2_h1",                    
    file == "chromsize/schMedA2_h2.fa.chromsize" ~ "schMedA2_h2",                                    
    file == "chromsize/schMedS3_h1.fa.chromsize" ~ "schMedS3_h1",                    
    file == "chromsize/schMedS3_h2.fa.chromsize" ~ "schMedS3_h2",     
    file == "chromsize/schNov1.fa.chromsize" ~ "schNov1",                                    
    file == "chromsize/schLug1.fa.chromsize" ~ "schLug1",                    
    file == "chromsize/schPol2.fa.chromsize" ~ "schPol2" 
  ))

res_lst <- list()

res_lst[[1]] <- parse_full_table("01_busco/schLug1/busco/run_metazoa_odb10/full_table.tsv", "schLug1")
res_lst[[2]] <- parse_full_table("01_busco/schPol2//busco/run_metazoa_odb10/full_table.tsv", "schPol2")
res_lst[[3]] <- parse_full_table("01_busco/schMedA2_h1/busco/run_metazoa_odb10/full_table.tsv", "schMedA2_h1")
res_lst[[4]] <- parse_full_table("01_busco/schMedA2_h2/busco/run_metazoa_odb10/full_table.tsv", "schMedA2_h2")
res_lst[[5]] <- parse_full_table("01_busco/schMedS3_h1/busco/run_metazoa_odb10/full_table.tsv", "schMedS3_h1")
res_lst[[6]] <- parse_full_table("01_busco/schMedS3_h2/busco/run_metazoa_odb10/full_table.tsv", "schMedS3_h2")
res_lst[[7]] <-  parse_full_table("01_busco/clonorchis_sinensis.PRJNA386618/busco/run_metazoa_odb10/full_table.tsv", "cloSin")
res_lst[[8]] <-  parse_full_table("01_busco/hymenolepis_microstoma.PRJEB124/busco/run_metazoa_odb10/full_table.tsv", "hymMic")
res_lst[[9]] <-  parse_full_table("01_busco/schistosoma_mansoni.PRJEA36577/busco/run_metazoa_odb10/full_table.tsv", "schMan")
res_lst[[10]] <- parse_full_table("01_busco/taenia_multiceps.PRJNA307624/busco/run_metazoa_odb10/full_table.tsv", "taeMul")
res_lst[[11]] <- parse_full_table("01_busco/braLan3/busco/run_metazoa_odb10/full_table.tsv", "braLan3")
res_lst[[12]] <- parse_full_table("01_busco/schNov1/busco/run_metazoa_odb10/full_table.tsv", "schNov1")

# add chromsize
res_lst <- lapply(res_lst, left_join, all_chromsize, by = c('sp', 'sequence'))
# braLan3 amphioxus
# chromosome names in these assemblies vary quite a bit
# I am here interested in the large contigs that correspond to the chromosomes and I will discard the small contigs
# same goes for haplotigs. For the broad scale analysis we do here they are not needed.

# schPol2
res_lst[[2]] <- res_lst[[2]] %>% 
  filter(grepl("chr[1234]$", sequence))
# schMedA2_h1
res_lst[[3]] <- res_lst[[3]] %>% 
  filter(grepl("Chr", sequence)) %>%
  mutate(sequence = gsub("Chr([0-9])_h.", "chr\\1", sequence))
# schMedA2_h2
res_lst[[4]] <- res_lst[[4]] %>% 
  filter(grepl("Chr", sequence)) %>%
  mutate(sequence = gsub("Chr([0-9])_h.", "chr\\1", sequence))
# schMedS3_h1
res_lst[[5]] <- res_lst[[5]] %>% 
  filter(grepl("chr", sequence)) %>%
  mutate(sequence = gsub("chr([0-9])_h.", "chr\\1", sequence))
# schMedS3_h2
res_lst[[6]] <- res_lst[[6]] %>% 
  filter(grepl("chr", sequence)) %>%
  mutate(sequence = gsub("chr([0-9])_h.", "chr\\1", sequence))

# cloSin has 28 chromosomes.trematoda
# I will here take only the numbered ones
res_lst[[7]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst[[7]] <- res_lst[[7]] %>% 
  filter(sequence %in% 1:7) %>%
  mutate(sequence = paste0("chr", sequence))

# hymMic cestoda, six main chromosomes
res_lst[[8]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst[[8]] <- res_lst[[8]] %>% 
  mutate(sequence = gsub(pattern = 'HMN_0(.)_pilon', replacement = 'chr\\1', sequence))

# schMan trematode with six autosomes and a ZW sex chromosome
# Assembly contains haplotypes marked with "H
res_lst[[9]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst[[9]] <- res_lst[[9]] %>% 
  filter(! (grepl("H0", sequence) | grepl("U0", sequence))) %>%
  mutate(sequence = gsub(pattern = 'SM_V7_', replacement = 'chr', sequence))

# taeMul cestoda, not sure about number of chromosomes but there are only seven that are hit.
res_lst[[10]] %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n())

res_lst[[10]] <- res_lst[[10]] %>% 
  mutate(sequence = gsub(pattern = 'LG', replacement = 'chr', sequence))


# braLan amphioxous, chromosome scaffolds start with OV and seem to start with OV696686.1
res_lst[[11]] %>%
  group_by(sp, sequence) %>% 
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n()) %>%
  arrange(desc(number_of_hits)) 

res_lst[[11]] <- res_lst[[11]] %>%
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

## Basic summary stats

################################################################################

all %>%
  group_by(sp, sequence) %>%
  filter(! (status == "Missing")) %>%
  summarise(number_of_hits = n()) %>%
  View()


################################################################################

## only S3_h1

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

## version with schMan as the source of the color

################################################################################
col_source <- all |> 
  filter(sp == 'schMan') |> 
  arrange(sequence, start)

pl <- setNames(RColorBrewer::brewer.pal(n=8, 'Dark2'), c(paste0('chr', 1:7), 'chrZW'))

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


pdf('fig/busco_syn_flatworms_normalized_schMan_colors_grey_v2.pdf', width = 8, height = 10)

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
  geom_segment(aes(y = as.numeric(plotting_sp) - 0.2 - 1,
                   yend= as.numeric(plotting_sp) + 0.2 - 1,
                   xend=normalize_start)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        line = element_blank()) 

dev.off()




pdf('fig/busco_syn_flatworms_normalized_schMan_colors.pdf', width = 8, height = 10)

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
