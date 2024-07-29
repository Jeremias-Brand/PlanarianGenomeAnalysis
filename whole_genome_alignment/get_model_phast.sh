#!/bin/bash
# use the docker container for easy installation
# docker run -v $(pwd):/data -it  --rm quay.io/comparative-genomics-toolkit/cactus
MAF='schmidtea_Dez22_asex.maf'
REF=schMedS3h1
CHR=$1
REFGFF="2022-sch-comp/annotation/schMedS3_h1_ENCODE_hybrid_longest_coding.gff3.gz"
# trying with taffy
# the maf file is reference based 
# when we extract chromosomes we need to only look at the reference

# index the MAF first
# taffy index -i $MAF

mkdir -p taffy_phast
mkdir -p taffy_phast/chunks

# make the annotation file for the current chromsome
# MAKE SURE prepare_phast was run
# zgrep ${CHR} ${REFGFF} | sed "s/${CHR}/schMedS3h1.${CHR}/" > taffy_phast/${REF}_${CHR}.gff3
# # converst to gtf
# gffread taffy_phast/${REF}_${CHR}.gff3 -T  > taffy_phast/${REF}_${CHR}.gtf
# cut -f1,2 2022-sch-comp/cactus/schMedS3_h1.fa.masked.fai > taffy_phast/schMedS3_h1.chromsize



# Configuration
CHUNK_SIZE=10000000  # size of each chunk
CHROMSIZE_FILE=taffy_phast/schMedS3_h1.chromsize
MAX_JOBS=40  # Maximum number of parallel jobs

# Function to process a single chunk
process_chunk() {
    chrom=$1
    start_pos=$2
    end_pos=$3
    REF=$4
    MAF=$5

    OUTNAME="taffy_phast/chunks/${REF}.${chrom}_${start_pos}-${end_pos}"
    echo "Processing chunk ${chrom}:${start_pos}-${end_pos}"

    # Taffy operations
    taffy view --inputFile $MAF --region ${REF}.${chrom}:${start_pos}-${end_pos} > ${OUTNAME}.taf
    taffy norm --filterGapCausingDupes --halFile schmidtea_Dez22.hal --inputFile ${OUTNAME}.taf > ${OUTNAME}_norm.taf
    taffy view --maf --inputFile ${OUTNAME}_norm.taf --outputFile ${OUTNAME}_norm.maf

    # phylofit
    msa_view  ${OUTNAME}_norm.maf --in-format MAF  --4d --features taffy_phast/${REF}_${CHR}.gtf > ${OUTNAME}_norm_4d.ss 2> ${OUTNAME}_norm_4d.log
    # make the tuple size one for the model
    msa_view ${OUTNAME}_norm_4d.ss  --in-format SS --out-format SS --tuple-size 1 > ${OUTNAME}_norm_4d_t1.ss 2> ${OUTNAME}_norm_4d_t1.log

    # phylogenetic tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))"
    # nonconserved model
    phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --out-root ${OUTNAME}_norm_4d ${OUTNAME}_norm_4d_t1.ss 2> ${OUTNAME}_norm_4d_phylofit.log

    # get models from each codon
    msa_view ${OUTNAME}_norm.maf --features taffy_phast/${REF}_${CHR}.gtf \
    --catmap "NCATS = 3; CDS 1-3" --out-format SS --unordered-ss  > ${OUTNAME}_norm_codons.ss 2> ${OUTNAME}_norm_codons.log

    # conserved model with the first position
    phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --do-cats 1 \
            --out-root ${OUTNAME}_conserved-codon1 ${OUTNAME}_norm_codons.ss 2> ${OUTNAME}_norm_codons_phylofit.log

}

# Read chromosome sizes and process each in chunks
while IFS=$'\t' read -r chrom size; do
    # Only process the chromosome of interest
    if [[ $chrom == $CHR ]]; then
        end_pos=0
        job_count=0
        while [[ $end_pos -lt $size ]]; do
            start_pos=$end_pos
            ((end_pos=start_pos+CHUNK_SIZE))
            if [[ $end_pos -gt $size ]]; then
                end_pos=$size
            fi

            # Process in the background
            process_chunk $chrom $start_pos $end_pos $REF $MAF &
            ((job_count++))

            # Limit the number of parallel jobs
            if [[ $job_count -ge $MAX_JOBS ]]; then
                wait -n  # Wait for any job to finish
                ((job_count--))
            fi
        done

        # Wait for all jobs from this chromosome to finish
        wait
    fi
done < "$CHROMSIZE_FILE"

exit

msa_view --unordered-ss --out-format SS --aggregate schMedS3h1,schMedS3h2,schLug1,schNov1,schPol2 taffy_phast/chunks/schMedS3h1*4d_t1.ss | head


#!/bin/bash

# Configuration
CHUNK_SIZE=10000000  # size of each chunk
CHROMSIZE_FILE=taffy_phast/schMedS3_h1.chromsize

# Read chromosome sizes and process each in chunks
while IFS=$'\t' read -r chrom size; do
    # Only process the chromosome of interest
    if [[ $chrom == $CHR ]]; then
        end_pos=0
        while [[ $end_pos -lt $size ]]; do
            start_pos=$end_pos
            ((end_pos=start_pos+CHUNK_SIZE))
            if [[ $end_pos -gt $size ]]; then
                end_pos=$size
            fi

            OUTNAME="taffy_phast/chunks/${REF}.${chrom}_${start_pos}-${end_pos}"
            echo "Processing chunk from $start_pos to $end_pos"

            # Taffy operations
            taffy view --inputFile $MAF --region ${REF}.${chrom}:${start_pos}-${end_pos} > ${OUTNAME}.taf
            taffy norm --filterGapCausingDupes --halFile schmidtea_Dez22.hal --inputFile ${OUTNAME}.taf > ${OUTNAME}_norm.taf
            taffy view --maf --inputFile ${OUTNAME}_norm.taf --outputFile ${OUTNAME}_norm.maf

            # phylofit
            msa_view  ${OUTNAME}_norm.maf --in-format MAF  --4d --features taffy_phast/${REF}_${CHR}.gtf > ${OUTNAME}_norm_4d.ss
            # make the tuple size one for the model
            msa_view ${OUTNAME}_norm_4d.ss  --in-format SS --out-format SS --tuple-size 1 > ${OUTNAME}_norm_4d_t1.ss

            # phylogenetic tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))"
            # nonconserved model
            phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))"  --msa-format SS --out-root ${OUTNAME}_norm_4d ${OUTNAME}_norm_4d_t1.ss

            # get models from each codon
            msa_view ${OUTNAME}_norm.maf --features taffy_phast/${REF}_${CHR}.gtf \
            --catmap "NCATS = 3; CDS 1-3" --out-format SS --unordered-ss  > ${OUTNAME}_norm_codons.ss 

            # conserved model with the first position
            phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --do-cats 1 \
                    --out-root ${OUTNAME}_conserved-codon1 ${OUTNAME}_norm_codons.ss 

            # nonconserved model
            phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --out-root ${OUTNAME}_norm_catmap  ${OUTNAME}_norm_catmap.ss 

        done
    fi
done < "$CHROMSIZE_FILE"


exit












## code to define chunks to be extracted with taffy
## form of command
OUTNAME=taffy_phast/chunks/${REF}.${CHR}_0-10000000
taffy view --inputFile $MAF --region ${REF}.${CHR}:0-10000000 > ${OUTNAME}.taf


## we then process each of these chuncks
# normalize the maf file to reduce gaps
taffy norm --filterGapCausingDupes --halFile schmidtea_Dez22.hal  --inputFile ${OUTNAME}.taf > ${OUTNAME}_norm.taf
# convert to maf
taffy view --maf --inputFile ${OUTNAME}_norm.taf --outputFile ${OUTNAME}_norm.maf

# phast cons model estimation
# get the 4d model
msa_view  ${OUTNAME}_norm.maf --in-format MAF  --4d --features taffy_phast/${REF}_${CHR}.gtf > ${OUTNAME}_norm_4d.ss
# make the tuple size one for the model
msa_view ${OUTNAME}_norm_4d.ss  --in-format SS --out-format SS --tuple-size 1 > ${OUTNAME}_norm_4d_t1.ss

# phylogenetic tree
# "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))"
# nonconserved model
phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))"  --msa-format SS --out-root ${OUTNAME}_norm_4d ${OUTNAME}_norm_4d_t1.ss

# get models from each codon
# then get the first position
msa_view ${OUTNAME}_norm.maf --features taffy_phast/${REF}_${CHR}.gtf \
--catmap "NCATS = 3; CDS 1-3" --out-format SS --unordered-ss  > ${OUTNAME}_norm_codons.ss 

# conserved model
phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --do-cats 1 \
        --out-root ${OUTNAME}_conserved-codon1 ${OUTNAME}_norm_codons.ss 

# nonconserved model

phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --out-root ${OUTNAME}_norm_catmap  ${OUTNAME}_norm_catmap.ss 








msa_split $MAF --in-format MAF --refseq cactus/schMedS3_h1.fa.masked \
                --windows 1000000,0 --out-root CHUNKS/dd --out-format SS \
                --min-informative 1000 --between-blocks 5000 


## NOTES
taffy index -i ${OUTNAME}.taf
taffy stats -s -i ${OUTNAME}.taf
wc -l taffy_phast/schMedS3h1.chr4_h1.taf
55278593 taffy_phast/schMedS3h1.chr4_h1.taf