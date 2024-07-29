
# conda activate genespace
# export PATH=$PATH:/home/jbrand/tools/MCSanX
# R
library(GENESPACE)

runwd <- file.path('./genespace')
gpar <- init_genespace(
  wd = runwd, 
  path2mcscanx = "./tools/MCScanX",
  nCores = 144,
  orthofinderInBlk = TRUE)
gpar <- run_genespace(gsParam = gpar) 

gids <- c(
	"schLug1", 
	"schMedS3h1",
	"schMedS3h2",
	"schNov1",
	"schPol2"
)


parse_annotations(genespaceWd = runwd,
rawGenomeRepo = './genespace/rawGenomes',
genomeDirs = gids,
faString = 'pep',
gffString = 'gff3',
gffIdColumn = "ID",
gffStripText = "transcript:",
headerStripText = "",
headerEntryIndex = 1,
troubleShoot = TRUE)

out <- run_genespace(gsParam = gpar, overwrite = TRUE)