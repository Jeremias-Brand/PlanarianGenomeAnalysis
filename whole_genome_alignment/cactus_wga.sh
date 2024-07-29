# make whole genome alignment with cactus
docker pull quay.io/comparative-genomics-toolkit/cactus
docker run -v $(pwd):/data --rm quay.io/comparative-genomics-toolkit/cactus cactus jobStore_schmidtea \
cactus/cactus_config.txt schmidtea_Dez22.hal --logFile log__schmidtea_Dez22 &> output__schmidtea_Dez22


## Creating pairwise hals against the reference

```bash
DOCKER="docker run -v $(pwd):/data --rm quay.io/comparative-genomics-toolkit/cactus"

GENOMES=$( $DOCKER halStats --genomes schmidtea_Dez22.hal )
mkdir -p schmidtea_Dez22
for g in in $GENOMES;do
	mkdir -p schmidtea_Dez22/$g
	
$DOCKER halStats --genomes schmidtea_Dez22.hal
```
I move it into the base directory of the gbdb

## Phylop

```bash
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schMedS3h1 nonconserved_all_UCSC.mod schMedS3h1_phylop.wig &
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schMedS3h2 nonconserved_all_UCSC.mod schMedS3h2_phylop.wig &
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schMedA2h1 nonconserved_all_UCSC.mod schMedA2h1_phylop.wig &
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schMedA2h2 nonconserved_all_UCSC.mod schMedA2h2_phylop.wig &
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schMedS3h1 nonconserved_all_UCSC.mod schMedS3h1_phylop.wig &
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schLug1 nonconserved_all_UCSC.mod schLug1_phylop.wig &
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schNov1 nonconserved_all_UCSC.mod schNov1_phylop.wig &
halPhyloP --hdf5InMemory schmidtea_Dez22.hal schPol2 nonconserved_all_UCSC.mod schPol2_phylop.wig &
```

## Snps

```bash
halSnps --hdf5InMemory schmidtea_Dez22.hal schMedS3h1 schMedS3h2 --tsv schMedS3h1_schMedS3h2_snps.tsv &
halSnps --hdf5InMemory schmidtea_Dez22.hal schMedS3h1 schMedA2h1 --tsv schMedS3h1_schMedA2h1_snps.tsv &
halSnps --hdf5InMemory schmidtea_Dez22.hal schMedS3h1 schMedA2h2 --tsv schMedS3h1_schMedA2h2_snps.tsv & 
halSnps --hdf5InMemory schmidtea_Dez22.hal schMedA2h1 schMedA2h2 --tsv schMedA2h1_schMedA2h2_snps.tsv &
```

# Export cactus masked genomes

```bash
hal2fasta --hdf5InMemory schmidtea_Dez22.hal schMedS3h1 > schMedS3h1_cactus_masked.fa &
hal2fasta --hdf5InMemory schmidtea_Dez22.hal schMedS3h2 > schMedS3h2_cactus_masked.fa &
hal2fasta --hdf5InMemory schmidtea_Dez22.hal schMedA2h1 > schMedA2h1_cactus_masked.fa &
hal2fasta --hdf5InMemory schmidtea_Dez22.hal schMedA2h2 > schMedA2h2_cactus_masked.fa &
```



# Liftover

```bash
gff2bed < schMedS3_h1_ENCODE_hybrid_agat_only_hconf.gff3 > schMedS3_h1_only_hconf.bed

halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_only_hconf.bed schMedS3h1 liftover_schMedS3_h1_only_hconf_to_schMedS3h1.bed &
halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_only_hconf.bed schMedS3h2 liftover_schMedS3_h1_only_hconf_to_schMedS3h2.bed &
halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_only_hconf.bed schMedA2h1 liftover_schMedS3_h1_only_hconf_to_schMedA2h1.bed &
halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_only_hconf.bed schMedA2h2 liftover_schMedS3_h1_only_hconf_to_schMedA2h2.bed &

cp annotation/schMedS3_h1_ENCODE_hybrid_agat_hconf.gff3.gz .
gunzip schMedS3_h1_ENCODE_hybrid_agat_hconf.gff3.gz

gff2bed < schMedS3_h1_ENCODE_hybrid_agat_hconf.gff3 > schMedS3_h1_hconf.bed
halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_hconf.bed schMedS3h1 liftover_schMedS3_h1_hconf_to_schMedS3h1.bed &
halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_hconf.bed schMedS3h2 liftover_schMedS3_h1_hconf_to_schMedS3h2.bed &
halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_hconf.bed schMedA2h1 liftover_schMedS3_h1_hconf_to_schMedA2h1.bed &
halLiftover --hdf5InMemory --bedType 3 schmidtea_Dez22.hal schMedS3h1 schMedS3_h1_hconf.bed schMedA2h2 liftover_schMedS3_h1_hconf_to_schMedA2h2.bed &
```
