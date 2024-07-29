#!/bin/bash
# use the docker container for easy installation
# docker run -v $(pwd):/data -it  --rm quay.io/comparative-genomics-toolkit/cactus
MAF=$1
REF=$2
CHR=$3
MAF='schmidtea_Dez22_asex.maf'
REF=schMedS3h1
CHR=chr4_h1
mkdir -p tmp/
mafExtractor --maf $MAF --seq ${REF}.${CHR}  > tmp/${REF}.${CHR}.maf
# https://github.com/CshlSiepelLab/phast/issues/10
# split by chromosome
cut -f1,2 2022-sch-comp/cactus/schMedS3_h1.fa.masked.fai | \
awk -v maf="$MAF" -v ref="$REF" '{print "mafExtractor --maf " maf " --seq " ref "." $1 " --start 0 --stop " $2 " > tmp/" ref "-" $1 ".maf"}' | grep 'chr' > extract.cmd

bash extract.cmd

# split and conquer approach
# http://compgen.cshl.edu/phast/phastCons-HOWTO.html

# category for each codon
msa_view tmp/schMedS3h1-chr4_h1.maf --features tmp/schMedS3h1_chr4_h1.gff3 \
--catmap "NCATS = 3; CDS 1-3" --out-format SS > tmp/schMedS3h1_chr4_h1_catmap.ss 2> tmp/schMedS3h1_chr4_h1_catmap.log 

# 4-fold degenerate sites
msa_view  tmp/schMedS3h1_chr4_h1_catmap.ss --in-format SS --4d --features tmp/schMedS3h1_chr4_h1.gtf > tmp/schMedS3h1-chr4_h1_4d_catmap.ss
# make the tuple size one for the model
msa_view tmp/schMedS3h1-chr4_h1_4d_catmap.ss  --in-format SS --out-format SS --tuple-size 1 > tmp/schMedS3h1-chr4_h1_4d_catmap_t1.ss

# make a model out of it
phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --out-root tmp/schMedS3h1-chr4_h1_4d_catmap tmp/schMedS3h1-chr4_h1_4d_catmap_t1.ss

# use the model in phastcons
phastCons --target-coverage 0.25 --expected-length 12 --rho 0.4 \
        --most-conserved tmp/schMedS3h1-chr4_h1_most_conserved.bed --msa-format MAF tmp/schMedS3h1-chr4_h1.maf \
        tmp/schMedS3h1-chr4_h1_4d_catmap.mod > tmp/schMedS3h1-chr4_h1_phastcons.wig

# soft masked alignments are used. we might want to change that
msa_view tmp/schMedS3h1-chr4_h1.maf  --in-format MAF --out-format SS  > tmp/schMedS3h1-chr4_h1.ss

# phast cons with SS as input
phastCons --target-coverage 0.25 --expected-length 12 --rho 0.4 \
        --most-conserved tmp/schMedS3h1-chr4_h1_most_conserved_ss.bed --msa-format SS tmp/schMedS3h1-chr4_h1.ss \
        tmp/schMedS3h1-chr4_h1_4d_catmap.mod > tmp/schMedS3h1-chr4_h1_phastcons_ss.wig

# from an index file we can get the chr variables
# get nono-redundant annotation
zgrep 'chr4_h1' 2022-sch-comp/annotation/schMedS3_h1_ENCODE_hybrid_longest_coding.gff3.gz | sed 's/chr4_h1/schMedS3h1.chr4_h1/' > tmp/schMedS3h1_chr4_h1.gff3
# converst to gtf
gffread tmp/schMedS3h1_chr4_h1.gff3 -T  > tmp/schMedS3h1_chr4_h1.gtf

msa_view  tmp/schMedS3h1-chr4_h1.maf --in-format MAF --out-format SS --4d --features tmp/schMedS3h1_chr4_h1.gtf > tmp/schMedS3h1_chr4_h1_4d.ss


grep 'chr4_h1' 2022-sch-comp/schMedS3_h1_ENCODE_hybrid_agat_hconf.gff3 | sed 's/chr4_h1/schMedS3h1.chr4_h1/' > tmp/schMedS3h1_chr4_h1.gff3
gffread tmp/schMedS3h1_chr4_h1.gff3 -T  > tmp/schMedS3h1_chr4_h1.gtf

msa_view  tmp/schMedS3h1-chr4_h1.maf --in-format MAF --4d --features tmp/schMedS3h1_chr4_h1.gtf > tmp/schMedS3h1-chr4_h1_4d.ss

msa_view  tmp/schMedS3h1-chr4_h1.maf --in-format MAF --4d --features tmp/schMedS3h1_chr4_h1.gff3 > tmp/schMedS3h1-chr4_h1_4d.ss



msa_view tmp/schMedS3h1-chr4_h1_4d.ss --features tmp/schMedS3h1_chr4_h1.gff3 --4d --in-format SS --out-format SS --tuple-size 1 > tmp/schMedS3h1-chr4_h1_4d_t1.ss


--out-format SS
        --unordered-stats --tuple-size 3 --reverse-groups transcript_i

        --soft-masked
# msa_view tmp/schMedS3h1-chr4_h1_4d.ss  --in-format SS --out-format SS --tuple-size 1 > tmp/schMedS3h1-chr4_h1_4d_t1.ss

# http://compgen.cshl.edu/phast/phyloFit-tutorial.php#bcl

phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --out-root tmp/schMedS3h1-chr4_h1_4d tmp/schMedS3h1-chr4_h1_4d_t1.ss

phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --out-root tmp/schMedS3h1-chr4_h1_catmap tmp/schMedS3h1_chr4_h1_catmap.ss 

--most-conserved 


for i in {1..10};
do
msa_view chr${i}_3.sort.maf --in-format MAF --4d --features Zea_mays.chr$i.T01_cds.gtf >chr${i}_3.4d-codons.ss
msa_view chr${i}_3.4d-codons.ss --in-format SS --out-format SS --tuple-size 1 >chr${i}_3.4d-sites.ss
done
msa_view --aggregate maize,sorgum,millot chr1_3.4d-sites.ss chr2_3.4d-sites.ss chr3_3.4d-sites.ss chr4_3.4d-sites.ss chr5_3.4d-sites.ss chr6_3.4d-sites.ss 
chr7_3.4d-sites.ss chr8_3.4d-sites.ss chr9_3.4d-sites.ss chr10_3.4d-sites.ss >all.sites.ss





phyloFit -i MAF S3_A2_chr4.maf
phastCons --target-coverage 0.25 \
--expected-length 12 \
--rho 0.4 \
--msa-format MAF \
S3_A2_chr4.maf phyloFit.mod \
--seqname `cat $i | head -n 3  | tail -n 1 | tr -s ' ' | cut -f 2 -d ' ' | cut -d. -f2` \
--most-conserved mostCons/S3_A2_chr4.bed > wig/S3_A2_chr4.wig