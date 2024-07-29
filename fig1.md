
## Fig 1b

```bash
k=20
# collect these SRA datasets using fastq-dump
ids=(SRR1951777 SRR5408395 SRR954965 SRR959588 SRR959589)

for i in "${ids[@]}"; do
    echo "trimgalore --fastqc --gzip --paired $i_1.fastq.gz $i_2.fastq.gz --output_dir .";
done | parallel -j 20

for i in *.fq.gz;do echo "meryl count k=$k threads=10 memory=100 output ${i}_trim.meryl ${i}_trim.fq.gz > log__meryl_${i} 2>&1 "; done >> meryl1.cmd
for i in $( ls *_trim.fq.gz | cut -d_ -f1 | sort | uniq );do echo "meryl union-sum output ${i}.meryl ${i}_trim_*.meryl"; done >> meryl2.cmd
for a in *.fa; do for i in $( ls | grep meryl | grep -v fq );do echo "merqury.sh $i $a ${i%.meryl}_${a%.fa}" ; done & done >> meryl3.cmd

parallel -a meryl1.cmd -j 30 --joblog meryl1.log
bash meryl2.cmd
bash meryl3.cmd
```

Run `scripts/fig1_meryl.R`


## Fig 1c,1d


install [stained glass](https://github.com/mrvollger/StainedGlass)

Make alignment of genomes using minimap2 like described below.

Parse teh 

```R
# read the gap file

gaps <- read.table('dotplots/schMedS3_h1__schMedS2_guo_100kb_gaps_in_schMedS3_h1.tsv', 
                   header = TRUE)

################################################################################

## Make Bed files for each

################################################################################

################################################################################
## do each region individually


library(dplyr)

write_bed_and_fasta_commands <- function(df, genome_path, fa_path, config_path, bed_output_path) {
  # Create/open the file commands.txt and empty its contents if it already exists
  fileConn<-file("commands.txt")
  close(fileConn)
  
  for(i in 1:nrow(df)) {
    name <- paste0(df$seqnames[i], '_', df$start[i], '_', df$end[i])
    bed_file_path <- paste0(bed_output_path, name, ".bed")
    bed_content <- paste(df$seqnames[i], df$start[i], df$end[i], name, df$strand[i], sep="\t")
    write(bed_content, file=bed_file_path)
    
    command <- paste0("cat ", genome_path, " | seqkit subseq --bed ", bed_file_path, " > ", fa_path, name, ".fa")
    write(command, file="stained_glass/commands.cmd", append = TRUE) # append the command to the file
    
    # config files
    config_content <- paste(
      "sample: ", name, "\n",
      "fasta: ../", fa_path, name, ".fa\n",
      "window: 2000\n",
      "nbatch: 100\n",
      "alnthreads: 100\n",
      "mm_f: 10000\n",
      "tempdir: temp\n",
      sep=""
    )
    write(config_content, file=paste0(config_path, name, "_config.yaml"))
  }
}
genome_path <- "latest/schMedS3_h1.fa"
fa_path <- "stained_glass/fa/"
config_path <- "stained_glass/config/"
bed_output_path <- "stained_glass/bed/"

write_bed_and_fasta_commands(gaps, genome_path = genome_path,
                             fa_path = fa_path, bed_output_path = bed_output_path,
                             config_path = config_path)

```



```bash
#!/bin/bash

NAME=$1
GLASS=$2
CONFIGDIR="stained_glass/config/"
OUTDIR="stained_glass/results/"

cp ${CONFIGDIR}/${NAME}_config.yaml ${GLASS}/config/config.yaml
cd $GLASS
snakemake --use-conda --cores 12 && snakemake --use-conda --cores 12 make_figures && cd .. && mv ${GLASS}/results/${NAME}* $OUTDIR
```


## Fig 1e

#### Generate whole-genome alignments using minimap2

Subset to only chromosomes.

```bash
seqkit grep -r -p '"chr"' lib/schMedS3_h1.fa > lib/schMedS3_h1_chr.fa 
seqkit grep -r -p '"chr"' lib/schMedS3_h2.fa > lib/schMedS3_h2_chr.fa 
```
Align.

```bash
minimap2 -L -ax asm5 -t 10 --eqx lib/schMedS3_h1_chr.fa lib/schMedS3_h2_chr.fa | samtools sort -O BAM - > 08_wga/S3h1_S3h2_chr.bam 2> 08_wga/log__mm2_S3h1_S3h2_chr
```

### Synteny detection with syri

```bash
conda activate syri
export PATH=/projects/jere/dev/syri/syri/bin/:$PATH
syri -c S3h1_S3h2_chr.paf -r schMedS3_h1_chr.fa -q schMedS3_h2_chr.fa -F P --prefix S3h1_S3h2_chr
```

Filtering out smaller variants

```bash
bcftools view -i'TYPE!="snp"' 08_wga/S3h1_S3h2_chrsyri.vcf > 14_syri_filter/S3h1_S3h2_chrsyri_nosnp.vcf
bcftools view -i'TYPE!="indels"' 14_syri_filter/S3h1_S3h2_chrsyri_nosnp.vcf > 14_syri_filter/S3h1_S3h2_chrsyri_nosnp_noindel.vcf
cp 08_wga/S3h1_S3h2_chrsyri.vcf 14_syri_filter/

for f in $( ls 14_syri_filter/*.vcf  );do bgzip $f; done
for f in $( ls 14_syri_filter/*.vcf.gz  );do tabix -p vcf $f; done
bgzip 14_syri_filter/*.vcf
tabix -p vcf 14_syri_filter/*.vcf.gz
```

Plot the syri results.
`./scripts/fig1_syri.R`

## Fig 1g, 1j

```R
library(tidyverse)
library(ggpubr)
library(data.table)

tx_subset <- c( "dd_Smed_v6.PCF","dd_Smes_v1.PCF",
		"schMedS3_BH_all", "schMedS3_BH_hconf",
		"schMedS3_h1_all", "schMedS3_h1_hconf", "schMedS3_h2_all", "schMedS3_h2_hconf",
		"smes_v2_hconf_SMESG", "smes_v2_repeatfilt_SMESG", "Sm_Oxford_v1" )

### ############################################################################
### BUSCOs
### ############################################################################

buscos <- read_tsv("scripts_luca/BUSCO.stats.log", col_names=NA) %>%
		dplyr::select(X1, X3, X4) %>%
		rename(Transcriptome=X1, 
			N=X3, BUSCO=X4) %>%
		filter(Transcriptome %in% tx_subset)	

busco_ranking <- buscos %>% 
		filter(BUSCO!="Missing") %>% 
		group_by(Transcriptome) %>%
		summarise(N=sum(N)) %>%
		arrange(desc(N)) %>% 
		pull(Transcriptome)

buscos <- mutate(buscos, Transcriptome = factor(Transcriptome, levels=busco_ranking))
buscos |> 
#filter(buscos, BUSCO != "Missing") %>%
	ggbarplot("Transcriptome", "N", 
		fill="BUSCO", color="BUSCO", 
		label= TRUE, lab.col = "black", 
		lab.pos = "in", lab.size= 2.75,
		lab.vjust = 1, 
		palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
		theme(axis.text.x = element_text(angle = 90, 
				vjust = 0.5, hjust=1)) +
		ggtitle("BUSCO Completeness") +
		ylab("BUSCO Hits") +
		theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
			plot.margin = margin(1,1,1.5,1.2, "cm"))
		
ggsave("BUSCO.pdf", width=4, height=6)


# Genomic busco

parse_full_table <- function(tbl) {
  
  assembly_name <- gsub('/.+','', 
                        gsub('01_busco/','', tbl))
  
  h <- c("busco_id", "status", "sequence", "start", "end", "strand", "score", "length")
  busco <- fread(tbl, fill = TRUE, sep = "\t", skip = 3, header = FALSE)
  names(busco) <- h
  busco$assembly <- assembly_name
  return(busco)
}


tbls_files <- list.files('.', pattern = 'full_table', recursive = TRUE)
tbls <- lapply(tbls_files, parse_full_table)

df <- do.call(rbind, tbls) 
# the duplicated genes are present several times
# we need to filter them out
sum <- df |> 
  select(busco_id, status, assembly) |> 
  unique() |> 
  group_by(assembly, status) |> 
  summarise(N= dplyr::n(),
            percentage = 100*(dplyr::n()/954)) |> 
  rename(BUSCO = status)


sel = c('dd_Smes_g4',
        'schMedS2_guo',
        'schMedS3_h1',
        'schMedS3_h2')

okabe_ito_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sum|> 
  filter(assembly %in% sel) |> 
  #filter(buscos, BUSCO != "Missing") %>%
  ggbarplot("assembly", "N", 
            fill="BUSCO", color="BUSCO", 
            label= T, lab.col = "black", 
            lab.pos = "in", lab.size= 3.75,
            palette = c("#56B4E9",  "#009E73", "#F0E442",  "#D55E00")) +
  scale_x_discrete(labels = 
                     c('dd_Smes_g4',
                   'schMedS2',
                   'schMedS3h1',
                   'schMedS3h2') )+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1)) +
  ggtitle("BUSCO Completeness") +
  ylab("BUSCO Hits") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        plot.margin = margin(1,1,1.5,1.2, "cm"),
        legend.position = 'right',
        axis.title.x = element_blank()) 

ggsave("genome_BUSCO_smed.pdf", width = 120, height = 180, units = 'mm')


sum|> 
  filter(! assembly %in% c('dd_Smes_g4',
                           'schMedS3_h1h2',
                           'schMedS3_h1',
                           'schMedS3_h2',
                         'schMedA2_h1',
                         'schMedA2_h2')) |> 
  mutate(assembly = factor(assembly, levels = c('schPol2','schNov1', 'schLug1'))) |> 
  #filter(buscos, BUSCO != "Missing") %>%
  ggbarplot("assembly", "N", 
            fill="BUSCO", color="BUSCO", 
            label= TRUE, lab.col = "black", 
            lab.pos = "in", lab.size= 2.75,
            palette = c("#56B4E9",  "#00AFBB", "#E7B800", "#FC4E07")) +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1)) +
  ylab("BUSCO Hits") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        plot.margin = margin(1,1,1.5,1.2, "cm"))


ggsave("genome_BUSCO_other_schmidtea.pdf", width=4, height=6)


################################################################################

## reading the 96 genes

################################################################################


library(readxl)
library(dplyr)
library(tidyr)
library(forcats)
cand <- read.csv('fig1/candidates_96_plotting.csv') |> 
  na.omit()

p_df <- cand |> 
  select(- short_name, - Intron.Retention) |> 
  pivot_longer(-name, names_to = 'type', values_to = 'count') |> 
  mutate(type = factor(type, levels = c("Perfect", "Frame.Shift", "Truncation", "Chimera", "Fragmentation", "Missing")))


order <- p_df |> 
  filter(type == 'Perfect') |> 
  arrange(desc(count)) |> 
  pull(name)

p_df <- p_df |> 
  mutate(name = factor(name, levels = order))

p_df |> 
ggbarplot("name", "count", 
          fill="type", color="type", 
          label= TRUE, lab.col = "black", 
          lab.pos = "in", lab.size= 2.75,
          palette = c( '#00AFBB', '#4E60F5', "#E7B800", "#FC4E07", "#FE9405", 'grey')) +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 96)) +
  ggtitle("Annotation assessment of 96 Genes") +
  ylab("Number of Genes") +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        plot.margin = margin(1,1,1.5,1.2, "cm"),
        axis.text = element_text(size = 11))


ggsave("plot_96_genes.pdf", width=3, height=6)

```























