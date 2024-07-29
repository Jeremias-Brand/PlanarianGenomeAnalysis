#!/bin/bash
set -euo pipefail
DIRNAME=$1
GENOME=$2

OUTDIR=genespace/rawGenomes/${DIRNAME}/
mkdir -p ${OUTDIR}

cp latest_annotation/subsets/${GENOME}_ENCODE_hybrid_longest_coding_sorted.gff3 ${OUTDIR}
cp latest/pep/${GENOME}_longest_coding.pep ${OUTDIR}
