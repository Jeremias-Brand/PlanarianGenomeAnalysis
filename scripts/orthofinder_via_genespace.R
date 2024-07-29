library(GENESPACE)
runwd <- file.path('./genespace_parasites_metazoa')
gids <- c(
	"NemVec",
	"BraLan",
	"schLug1", 
	"schMedS3h1",
	"schMedS3h2",
	"schNov1",
	"schPol2",
	"cloSin",
	"echMul",
	"hymMic",
	"schMan",
	"taeMul"
)

gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  ploidy = rep(1, length(gids)),
  wd = runwd, 
  nCores = 144,
  nGaps = 5,
  blkSize = 2,
  blkRadius = 25,
  synBuff = 100,
  arrayJump = 50,
  orthofinderInBlk = TRUE,
  gffString = "gff3", 
  pepString = "pep",
  path2orthofinder = "orthofinder", 
  path2mcscanx = "~/tools/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

out <- run_genespace(gsParam = gpar, overwrite = TRUE)