################################################################################

## replotting odp reciprocal hits

################################################################################

library(tidyverse)
library(data.table)
library(patchwork)

process_data <- function(data, y_pos_column, y_scaf_column, sort_y = FALSE) {
  if (sort_y) {
    ordered_y_scaffolds <- data %>%
      group_by(!!sym(y_scaf_column)) %>%
      summarize(max_val = max(!!sym(y_pos_column), na.rm = TRUE)) %>%
      arrange(-max_val) %>%
      pull(!!sym(y_scaf_column))
    data[[y_scaf_column]] <- factor(data[[y_scaf_column]], levels = ordered_y_scaffolds)
    
    # filter out the small chromosomes
    fi = 5e6
    sel <- data %>%
      group_by(!!sym(y_scaf_column)) %>%
      summarize(max_val = max(!!sym(y_pos_column), na.rm = TRUE)) %>%
      arrange(-max_val) %>%
      filter(max_val > fi) %>%
      pull(!!sym(y_scaf_column))
    
    data <- data %>% 
      filter(!!sym(y_scaf_column) %in% sel)
  }
  
  return(data)
}

plot_data <- function(filepath, y_pos_column, y_scaf_column, sp, sort_y = FALSE) {
  
  alg <- fread(filepath)
  
  if (y_scaf_column == 'schMan_scaf') {
    alg[[y_scaf_column]] = gsub('SM_V9_', '', alg[[y_scaf_column]])
  }
  
  if (sort_y) {
    ordered_y_scaffolds <- alg %>%
      group_by(!!sym(y_scaf_column)) %>%
      summarize(max_val = max(!!sym(y_pos_column), na.rm = TRUE)) %>%
      arrange(-max_val) %>%
      pull(!!sym(y_scaf_column))
    alg[[y_scaf_column]] <- factor(alg[[y_scaf_column]], levels = ordered_y_scaffolds)
    
    ################################################################################
    ## filter out the small chromosomes
    fi = 5e6
    sel <- alg %>%
      group_by(!!sym(y_scaf_column)) %>%
      summarize(max_val = max(!!sym(y_pos_column), na.rm = TRUE)) %>%
      arrange(-max_val) %>%
      filter(max_val > fi) |> 
      pull(!!sym(y_scaf_column))
    
    alg <- alg |> 
      filter(!!sym(y_scaf_column) %in% sel)
  }
  
  p <- alg |> 
    mutate(al = ifelse(whole_FET > 0.05, 0.4, 1)) |> 
    ggplot(aes(y = get(y_pos_column), x = BCnS_LGs_pos)) +
    geom_point(aes(fill = color, alpha = al), shape = 21, stroke = NA, size = 0.8) +
    facet_grid(get(y_scaf_column) ~ BCnS_LGs_scaf, 
               scales = 'free', space = 'free') +
    xlab('MALG') +
    ylab(sp) +
    scale_x_continuous(expand = c(-0.02, 0)) +  # Adjust x-axis expand values
    scale_y_continuous(expand = c(0, -0.02))  + 
    scale_fill_identity() + 
    scale_alpha_identity() +
    theme_minimal() +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.1),
      axis.text = element_blank(),
      strip.text.x = element_text(angle = 90, hjust = 1, size = 5),
      strip.text.y = element_text(angle = 0, hjust = 1, size = 5)
    )
  
  return(p)
}


plot_pairwise <- function(filepath, y_pos_column, y_scaf_column, x_pos_column, x_scaf_column, y_sp, x_sp, sort_y = FALSE, FET = FALSE) {
  
  alg <- fread(filepath)
  
  if (x_scaf_column == 'schMan_scaf') {
    alg[[x_scaf_column]] = gsub('SM_V9_', '', alg[[x_scaf_column]])
  }
  
  if (y_scaf_column == 'schMan_scaf') {
    alg[[y_scaf_column]] = gsub('SM_V9_', '', alg[[y_scaf_column]])
  }
  
  if (sort_y) {
    ordered_y_scaffolds <- alg %>%
      group_by(!!sym(y_scaf_column)) %>%
      summarize(max_val = max(!!sym(y_pos_column), na.rm = TRUE)) %>%
      arrange(-max_val) %>%
      pull(!!sym(y_scaf_column))
    alg[[y_scaf_column]] <- factor(alg[[y_scaf_column]], levels = ordered_y_scaffolds)
    
    ################################################################################
    ## filter out the small chromosomes
    fi = 5e6
    sel <- alg %>%
      group_by(!!sym(y_scaf_column)) %>%
      summarize(max_val = max(!!sym(y_pos_column), na.rm = TRUE)) %>%
      arrange(-max_val) %>%
      filter(max_val > fi) |> 
      pull(!!sym(y_scaf_column))
    
    alg <- alg |> 
      filter(!!sym(y_scaf_column) %in% sel)
  }
  
  
  ################################################################################
  ## switch for viz of FET significance
  if (FET == FALSE) {
    p_dat <-  alg |> 
      mutate(al = 1) 
  } else {
    p_dat <-  alg |> 
      mutate(al = ifelse(whole_FET > 0.05, 0.4, 1))
  }
  p <- p_dat |> 
    ggplot(aes(y = get(y_pos_column), x = get(x_pos_column))) +
    geom_point(aes(alpha = al), fill = "grey25", shape = 21, stroke = NA, size = 0.8) +
    facet_grid(get(y_scaf_column) ~ get(x_scaf_column), 
               scales = 'free', space = 'free') +
    xlab(x_sp) +
    ylab(y_sp) +
    scale_x_continuous(expand = c(-0.02, 0)) +  # Adjust x-axis expand values
    scale_y_continuous(expand = c(0, -0.02))  + 
    scale_alpha_identity() +
    theme_minimal() +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.1),
      axis.text = element_blank(),
      strip.text.x = element_text(angle = 90, hjust = 1, size = 5),
      strip.text.y = element_text(angle = 0, hjust = 1, size = 5)
    )
  
  return(p)
}




# To plot for the first dataset:
p_a1 <- plot_data('00_odp/step2-figures/ALG-species_plots/BCnS_LGs_schMedS3h1_xy_reciprocal_best_hits.plotted.rbh',
                  'schMedS3h1_pos', 'schMedS3h1_scaf', 'S. mediterranea')

# To plot for the second dataset:
p_a2 <- plot_data('00_odp/step2-figures/ALG-species_plots/BCnS_LGs_schMan_xy_reciprocal_best_hits.plotted.rbh',
                  'schMan_pos', 'schMan_scaf', 'S. mansoni')


# To plot for the second dataset:
p_a3 <- plot_data('00_odp/step2-figures/ALG-species_plots/BCnS_LGs_Machtxv2_xy_reciprocal_best_hits.plotted.rbh',
                  'Machtxv2_pos', 'Machtxv2_scaf', 'M. hystrix', sort_y = TRUE)


r1 <- plot_grid(p_a1, p_a2, p_a3, nrow = 1, align = 'h', labels = c('D', 'E', 'F'))

################################################################################

## make the supplementary plots

################################################################################

p_a1 <- plot_data('00_odp/step2-figures/ALG-species_plots/BCnS_LGs_hymMic_xy_reciprocal_best_hits.plotted.rbhh',
                  'hymMic_pos', 'hymMic_scaf', 'S. mediterranea')



################################################################################

## replot pairwise synteny without color

################################################################################




p_s1 <- plot_pairwise('00_odp/step2-figures/synteny_nocolor/schMan_schMedS3h1_xy_reciprocal_best_hits.plotted.rbh',
              x_pos_column = 'schMedS3h1_pos', x_scaf_column = 'schMedS3h1_scaf', x_sp = 'S. mediterranea',
              y_pos_column = 'schMan_pos', y_scaf_column = 'schMan_scaf', y_sp = 'S. mansoni')

p_s2 <- plot_pairwise('00_odp/step2-figures/synteny_nocolor/Machtxv2_schMedS3h1_xy_reciprocal_best_hits.plotted.rbh',
              x_pos_column = 'schMedS3h1_pos', x_scaf_column = 'schMedS3h1_scaf', x_sp = 'S. mediterranea',
              y_pos_column = 'Machtxv2_pos', y_scaf_column = 'Machtxv2_scaf', y_sp = 'M. hystrix',
              sort_y = TRUE)

p_s3 <- plot_pairwise('00_odp/step2-figures/synteny_nocolor/Machtxv2_schMan_xy_reciprocal_best_hits.plotted.rbh',
              x_pos_column = 'schMan_pos', x_scaf_column = 'schMan_scaf', x_sp = 'S. mansoni',
              y_pos_column = 'Machtxv2_pos', y_scaf_column = 'Machtxv2_scaf', y_sp = 'M. hystrix',
              sort_y = TRUE)


r2 <- plot_grid(p_s1, p_s2, p_s3, nrow = 1, align = 'v', labels = c('G', 'H', 'I'))


png('test.png', width = 9, height = 6, units = 'in', res = 300)

plot_grid(r1, r2, 
          ncol = 1)

dev.off()

  
################################################################################

## using the patchwork package to combine the plots

################################################################################

theme_patchwork <- theme(plot.tag.position = c(0.2, 1.005))

a_plots <- {p_a2 + labs(tag = 'D')} + {p_a1 + labs(tag = 'E')} + {p_a3 + labs(tag = 'F')}

s_plots <- {p_s1 + labs(tag = 'G')} + {p_s2 + labs(tag = 'H')} + {p_s3 + labs(tag = 'I')}

png('test.png', width = 9, height = 6, units = 'in', res = 300)

a_plots / s_plots + theme_patchwork

dev.off()


  
pdf('MALG.pdf', width = 9, height = 6)

a_plots / s_plots + theme_patchwork

dev.off()

  

################################################################################

## Replot the other species

################################################################################

################################################################################
## MALG

RBHs <- c(
  list.files('00_odp/step2-figures/ALG-species_plots/', pattern = '*rbh', full.names = TRUE),
  '00_odp/step2-figures_old/ALG-species_plots/BCnS_LGs_cloSin_xy_reciprocal_best_hits.plotted.rbh'
)
path = '00_odp/step2-figures/ALG-species_plots/BCnS_LGs_schMan_xy_reciprocal_best_hits.plotted.rbh'
plot_wrapper <- function(path) {
  name = gsub('.+/BCnS_LGs_([^_]+)_xy_reciprocal_best_hits.plotted.rbh', '\\1', path)
  return(
    plot_data(
    path,
    paste0(name, '_pos'),
    paste0(name, '_scaf'),
    name
  )
  )
}


pdf('fig/odp/MALG_supplement.pdf', width = 6, height = 6)
lapply(RBHs, plot_wrapper)
dev.off()

system('gs -sDEVICE=png16m -o fig/odp/MALG_supp%03d.png -r600 fig/odp/MALG_supplement.pdf')

################################################################################
## only plot selection
selection <- data.frame(file = RBHs) |> 
  filter(
    grepl(paste(parasites, collapse='|'), file) |
      grepl(paste(schmidtea, collapse='|'), file) |
      grepl(paste(mac, collapse='|'), file) &
      !grepl('Maccliv2|schMedS3h2', file)) |> 
  pull(file)

sel_plots <- lapply(selection[c(5,7,6,3, 9,4, 1,8,2)], plot_wrapper)
x_plots <- {sel_plots[[1]] + labs(tag = 'A')} + {sel_plots[[2]] + labs(tag = 'B')} + {sel_plots[[3]] + labs(tag = 'C')}

y_plots <- {sel_plots[[4]] + labs(tag = 'D')} + {sel_plots[[5]] + labs(tag = 'E')} + {sel_plots[[6]] + labs(tag = 'F')}

z_plots <- {sel_plots[[7]] + labs(tag = 'G')} + {sel_plots[[8]] + labs(tag = 'H')} + {sel_plots[[9]] + labs(tag = 'I')}
final_plot <- x_plots / y_plots / z_plots + theme_patchwork

ggsave(filename = 'fig/odp/MALG_selection.pdf', plot = final_plot, width = 18, height = 18)

system('gs -sDEVICE=png16m -o fig/odp/MALG_selection.png -r600 fig/odp/MALG_selection.pdf')

################################################################################
## pairwise

cloSin <- list.files('00_odp/step2-figures_old/synteny_nocolor', pattern = '*cloSin*', full.names = TRUE)

RBHs <- c(
  list.files('00_odp/step2-figures/synteny_nocolor', pattern = '*rbh', full.names = TRUE),
  cloSin[grepl('.rbh', cloSin)]
)

pair_wrapper <- function(path) {
  x_name = gsub('.+/([^_]+)_([^_]+)_xy_reciprocal_best_hits.plotted.rbh', '\\1', path)
  y_name = gsub('.+/([^_]+)_([^_]+)_xy_reciprocal_best_hits.plotted.rbh', '\\2', path)
  return(
    plot_pairwise(
      path,
      x_pos_column = paste0(x_name,'_pos'),
      x_scaf_column = paste0(x_name,'_scaf'),
      x_sp = x_name,
      y_pos_column = paste0(y_name,'_pos'),
      y_scaf_column = paste0(y_name,'_scaf'),
      y_sp = y_name,)
  )
}

plot_list <- lapply(RBHs, pair_wrapper)
names(plot_list) <- RBHs

################################################################################
## composite plots
parasites <- c('cloSin', 'schMan', 'hymMic', 'taeMul')
schmidtea <- c('schMedS3h1', 'schPol2', 'schNov1', 'schLug1')
mac <- c('Machtx')

################################################################################
## parasites only
p_only <- data.frame(file = RBHs) |> 
  filter(
    grepl(paste(parasites, collapse='|'), file) &
    !grepl(paste(schmidtea, collapse='|'), file) &
    !grepl(paste(mac, collapse='|'), file) &
    !grepl('Maccliv2|schMedS3h2', file)) |> 
  pull(file)

pp_only <- lapply(p_only, pair_wrapper)
pdf('fig/odp/MALG_pairwise_parasites.pdf', width = 18, height = 12)

a_plots <- {pp_only[[1]] + labs(tag = 'A')} + {pp_only[[2]] + labs(tag = 'B')} + {pp_only[[3]] + labs(tag = 'C')}

s_plots <- {pp_only[[4]] + labs(tag = 'D')} + {pp_only[[5]] + labs(tag = 'E')} + {pp_only[[6]] + labs(tag = 'F')}

a_plots / s_plots + theme_patchwork
dev.off()

system('gs -sDEVICE=png16m -o fig/odp/MALG_pairwise_parasites.png -r600 fig/odp/MALG_pairwise_parasites.pdf')


s_only <- data.frame(file = RBHs) |> 
  filter(
    !grepl(paste(parasites, collapse='|'), file) &
      grepl(paste(schmidtea, collapse='|'), file) &
      !grepl(paste(mac, collapse='|'), file) &
      !grepl('Maccliv2|schMedS3h2', file)) |> 
  pull(file)

ss_only <- lapply(s_only, pair_wrapper)
pdf('fig/odp/MALG_pairwise_schmidtea.pdf', width = 18, height = 12)

a_plots <- {ss_only[[1]] + labs(tag = 'A')} + {ss_only[[2]] + labs(tag = 'B')} + {ss_only[[3]] + labs(tag = 'C')}

s_plots <- {ss_only[[4]] + labs(tag = 'D')} + {ss_only[[5]] + labs(tag = 'E')} + {ss_only[[6]] + labs(tag = 'F')}

a_plots / s_plots + theme_patchwork
dev.off()

system('gs -sDEVICE=png16m -o fig/odp/MALG_pairwise_schmidtea.png -r600 fig/odp/MALG_pairwise_schmidtea.pdf')


sp_only <- data.frame(file = RBHs) |> 
  filter(
      grepl(paste(parasites, collapse='|'), file) &
      grepl(paste(schmidtea, collapse='|'), file) &
      !grepl(paste(mac, collapse='|'), file) &
      !grepl('Maccliv2|schMedS3h2', file)) |> 
  pull(file)

ssp_only <- lapply(sp_only, pair_wrapper)

pdf('fig/odp/MALG_pairwise_sp1.pdf', width = 18, height = 18)

a_plots <- {ssp_only[[1]] + labs(tag = 'A')} + {ssp_only[[2]] + labs(tag = 'B')} + {ssp_only[[3]] + labs(tag = 'C')}

s_plots <- {ssp_only[[4]] + labs(tag = 'D')} + {ssp_only[[5]] + labs(tag = 'E')} + {ssp_only[[6]] + labs(tag = 'F')}

x_plots <- {ssp_only[[7]] + labs(tag = 'G')} + {ssp_only[[8]] + labs(tag = 'H')} + {ssp_only[[9]] + labs(tag = 'I')}
a_plots / s_plots / x_plots + theme_patchwork
dev.off()

system('gs -sDEVICE=png16m -o fig/odp/MALG_pairwise_sp1.png -r600 fig/odp/MALG_pairwise_sp1.pdf')

pdf('fig/odp/MALG_pairwise_sp2.pdf', width = 18, height = 18)

a_plots <- {ssp_only[[10]] + labs(tag = 'A')} + {ssp_only[[11]] + labs(tag = 'B')} + {ssp_only[[12]] + labs(tag = 'C')}

s_plots <- {ssp_only[[13]] + labs(tag = 'D')} + {ssp_only[[14]] + labs(tag = 'E')} + {ssp_only[[15]] + labs(tag = 'F')}

x_plots <- {ssp_only[[16]] + labs(tag = 'G')} + plot_spacer() + plot_spacer()
a_plots / s_plots / x_plots + theme_patchwork
dev.off()

system('gs -sDEVICE=png16m -o fig/odp/MALG_pairwise_sp2.png -r600 fig/odp/MALG_pairwise_sp2.pdf')


################################################################################
## machtx to everybody

m_only <- data.frame(file = RBHs) |> 
  filter(
      grepl(paste(mac, collapse='|'), file) &
      !grepl('Maccliv2|schMedS3h2', file)) |> 
  pull(file)

mm_only <- lapply(m_only, pair_wrapper)



a_plots <- {mm_only[[1]] + labs(tag = 'A')} + {mm_only[[2]] + labs(tag = 'B')} + {mm_only[[3]] + labs(tag = 'C')}

s_plots <- {mm_only[[4]] + labs(tag = 'D')} + {mm_only[[5]] + labs(tag = 'E')} + {mm_only[[6]] + labs(tag = 'F')}

x_plots <- {mm_only[[7]] + labs(tag = 'G')} + {mm_only[[8]] + labs(tag = 'H')} + plot_spacer()

final_plot <- a_plots / s_plots / x_plots + theme_patchwork
ggsave(filename = 'fig/odp/MALG_pairwise_m.pdf', plot = final_plot, width = 18, height = 18)


system('gs -sDEVICE=png16m -o fig/odp/MALG_pairwise_m.png -r600 fig/odp/MALG_pairwise_m.pdf')



pdf('fig/odp/MALG_pairwise_supplement.pdf', width = 6, height = 6)
lapply(RBHs, pair_wrapper)
dev.off()

system('gs -sDEVICE=png16m -o fig/odp/MALG_pairwise_supp%03d.png -r600 fig/odp/MALG_pairwise_supplement.pdf')


