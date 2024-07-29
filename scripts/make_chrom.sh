#!/bin/bash
# agat can be installed from source, conda or use the docker: docker run --rm  -v `pwd`:/data quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0
set -euo pipefail
NAME=$1
DIR=$2
agat_convert_sp_gff2bed.pl -gff $DIR/${NAME}.gff3 -o $DIR/${NAME}.bed
awk 'BEGIN {OFS = "\t"}{print $4,$1,$6,$2,$3}' ${DIR}/${NAME}.bed > ${DIR}/${NAME}.chrom
