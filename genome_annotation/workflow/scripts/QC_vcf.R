# https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html
# https://cran.r-project.org/web/packages/PopGenome/vignettes/An_introduction_to_the_PopGenome_package.pdf
library(vcfR)
library(ggplot2)
library(ggthemes)
ref <- ape::read.dna("resources/dd_Smes_g4.fasta", format = "fasta")
chs <- names(dna)
# let's look at the first few
i=4
# add a space so it only takes exact matches
dna <- ref[ names(ref) == chs[i] ]
dna <- as.matrix(dna)
vcf <- read.vcfR(paste0("06_vcf/scaffolds/",chs[i],".vcf.gz"))
chrom <- create.chromR(chs[i], vcf=vcf, seq = dna)
plot(chrom)

# custom plots
vcf_gg <- extract_gt_tidy(vcf)

ggplot(vcf_gg, aes(x = gt_DP)) +
  geom_density(color="black", fill="white") +
  facet_wrap(~Indiv)


ggplot(vcf_gg, aes(x = gt_DP, fill = Indiv, color = Indiv)) +
  geom_density(alpha =0.2) +
  #facet_wrap(~Indiv) +
  theme_bw()

ggplot(vcf_gg, aes(x = Indiv, y = gt_DP, fill = Indiv, color = Indiv)) +
  geom_violin(alpha =0.2) +
  #facet_wrap(~Indiv) +
  theme_bw() + 
  scale_fill_colorblind() +
  scale_color_colorblind()

ggplot(vcf_gg, aes(x  = gt_DP, fill = Indiv, color = Indiv)) +
  geom_density(alpha =0.2) +
  theme_bw() + 
  scale_fill_colorblind() +
  scale_color_colorblind()


strwrap(vcf@meta[1:20])
queryMETA(vcf)
queryMETA(vcf, element = 'DP')
queryMETA(vcf, element = 'QD')
head(getFIX(vcf))


vcf@fix[,vcf@fix[,1] == chs[1]]
extract_gt_tidy(vcf)

vcf@gt[1:6, 1:5]

# first contig
chrom <- create.chromR("ch100", vcf=vcf, seq = dna)
plot(chrom)
chrom <- masker(chrom, min_QUAL = 1, min_DP = 10, max_DP = 700)
plot(chrom)
chrom_p <- proc.chromR(chrom, verbose=TRUE)
plot(chrom_p)

chrom_f <- proc.chromR(chrom_f, verbose=TRUE)
plot(chrom_f)
chromoqc(chrom, dp.alpha=20)
chromoqc(chrom, xlim=c(5e+05, 6e+05))



plot(chrom)

vcf@gt
vcf <- read.vcfR("06_vcf/all_raw.vcf.gz")
pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

library(vcfR)
vcf <- read.vcfR( vcf_file, verbose = FALSE )
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")


chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)


plot(chrom)


chrom_f <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)
plot(chrom_f)
chrom_f <- proc.chromR(chrom_f, verbose=TRUE)
plot(chrom_f)
chromoqc(chrom, dp.alpha=20)
chromoqc(chrom, xlim=c(5e+05, 6e+05))
