library(dplyr)
library(ggplot2)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Plot histogram of Khmer output")
# Add command line arguments
p <- add_argument(p, "--infile", help="report of khmer output", type="character")
p <- add_argument(p, "--outfile", help="output file", default="kmer_histo.pdf")
# Parse the command line arguments
argv <- parse_args(p)

df <- read.table(argv$infile)
names(df) <- c("genome", "kmer", "size")
genome = unique(df$genome)

pdf(argv$outfile)
ggplot(df, aes(x=as.factor(kmer), y = size/1e6)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ggtitle(paste0("Effective Genome Size of ", genome)) +
  xlab("khmer k-mer size") +
  ylab("Effective genome size (Mb)") +
  theme_classic()
dev.off()