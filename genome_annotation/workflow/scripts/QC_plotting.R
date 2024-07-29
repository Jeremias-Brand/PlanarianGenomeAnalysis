library(vcfR)
library(ggplot2)
library(ggthemes)


vcf <- read.vcfR("06_vcf/single/SRR3398520.vcf.gz")
vcf_gg <- extract_gt_tidy(vcf)
vcf_info_gg <- extract_info_tidy(vcf)


names(vcf_gg)
names(vcf_info_gg)

# plot Quality/depth QD
QD_filter = 4
ggplot(vcf_info_gg, aes(x = QD)) +
  geom_histogram(color="black", fill="steelblue") +
  geom_vline(xintercept = QD_filter, color = "red", size = 2, lty = 2) +
  theme_bw()


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
