#!/bin/bash
cd ..
DOCKER="docker run -v $(pwd):/data --rm -u $(id -u):$(id -g) quay.io/comparative-genomics-toolkit/cactus"
WGA=whole_genome_alignment/schmidtea_Dez22.hal
SOURCEDIR=conserved_elements/ATAC_peak_files/
SINKDIR=conserved_elements/schMedS3h1_based_full/
for i in schMedS3h2 schPol2 schNov1 schLug1; do 
$DOCKER halLiftover --bedType 4 $WGA schMedS3h1 ${SOURCEDIR}schMedS3h1_full_peak.tab \
$i ${SINKDIR}${i}_halLiftover_schMedS3h1_full_peak.bed;  done
