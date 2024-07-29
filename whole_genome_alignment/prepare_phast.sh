MAF='schmidtea_Dez22_asex.maf'
REF=schMedS3h1
CHR=$1
REFGFF="2022-sch-comp/annotation/schMedS3_h1_ENCODE_hybrid_longest_coding.gff3.gz"

zgrep ${CHR} ${REFGFF} | sed "s/${CHR}/schMedS3h1.${CHR}/" > taffy_phast/${REF}_${CHR}.gff3
# converst to gtf
gffread taffy_phast/${REF}_${CHR}.gff3 -T  > taffy_phast/${REF}_${CHR}.gtf
cut -f1,2 2022-sch-comp/cactus/schMedS3_h1.fa.masked.fai > taffy_phast/schMedS3_h1.chromsize