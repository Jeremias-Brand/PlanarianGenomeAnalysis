```bash

cp -r /DATA11900/projects/Comparative_Genomics_Schmidtea/rinklab_genomes/latest .
cp -r /DATA11900/projects/Comparative_Genomics_Schmidtea/rinklab_genomes/latest_annotation .

# we need pep and gff files in the correct folder structure
# the folders can't have special characters
./setup_genespace_directory.sh schLug1 schLug1
./setup_genespace_directory.sh schMedS3h1 schMedS3_h1
./setup_genespace_directory.sh schMedS3h2 schMedS3_h2
./setup_genespace_directory.sh schNov1 schNov1
./setup_genespace_directory.sh schPol2 schPol2

conda activate genespace
export PATH=$PATH:/home/jbrand/tools/MCSanX
```
Run interactively

```R
#devtools::install_github("jtlovell/GENESPACE@dev", upgrade = F)  
library(GENESPACE)

runwd <- file.path('$PWD/genespace')

gids <- c(
	"schLug1", 
	"schMedS3h1",
	"schMedS3h2",
	"schNov1",
	"schPol2")

parse_annotations(genespaceWd = runwd,
rawGenomeRepo = '$PWD/genespace/rawGenomes',
genomeDirs = gids,
faString = 'pep',
gffString = 'gff3',
gffIdColumn = "ID",
gffStripText = "ID=",
headerEntryIndex = 1,
troubleShoot = TRUE)

gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  ploidy = rep(1, length(gids)),
  wd = runwd, 
  nCores = 140,
  gffString = "gff3", 
  pepString = "pep",
  path2orthofinder = "orthofinder", 
  path2mcscanx = "~/tools/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))

out <- run_genespace(gsParam = gpar, overwrite = TRUE)
