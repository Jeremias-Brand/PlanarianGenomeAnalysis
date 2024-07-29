#!/bin/bash
# To setup follow code and instructions in ./setup/setup_odp.sh

wget https://zenodo.org/record/7861770/files/macrostomum_genomes_data_repo.tar.gz
tar xvf macrostomum_genomes_data_repo.tar.gz
mkdir -p genomes
mv macrostomum_genomes_data_repo/assemblies/Machtx_SR1_v2.fasta genomes/
mv macrostomum_genomes_data_repo/assemblies/Machtx_SR1_v2_rnm_rmnc_longest.gff3 genomes/
.scripts/make_chrom.sh Machtx_SR1_v2_rnm_rmnc_longest genomes
gffread -y genomes/Machtx_SR1_v2.pep -g genomes/Machtx_SR1_v2.fasta genomes/Machtx_SR1_v2_rnm_rmnc_longest.gff3

bash ./scripts/run_make_chrom.sh


sed -i 's/transcript://' cloSin.chrom
sed -i 's/transcript://' taeMul.chrom
sed -i 's/transcript://' hymMic.chrom
sed -i 's/transcript://' schMan.chrom


conda activate agat
NAME=schPol2
NAME=schLug1
NAME=schNov1
NAME=schMedS3_h1
NAME=schMedS3_h2
agat_convert_sp_gff2bed.pl -gff ann/${NAME}.gff3 -o ann/${NAME}.bed
awk 'BEGIN {OFS = "\t"}{print $4,$1,$6,$2,$3}' ann/${NAME}.bed > ann/${NAME}.chrom


