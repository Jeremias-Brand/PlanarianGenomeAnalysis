
msa_view --unordered-ss --out-format SS --aggregate schMedS3h1,schMedS3h2,schLug1,schNov1,schPol2 \
taffy_phast/chunks/schMedS3h1*4d_t1.ss > taffy_phast/schMedS3h1_all_4d.ss

msa_view --unordered-ss --out-format SS --aggregate schMedS3h1,schMedS3h2,schLug1,schNov1,schPol2 \
taffy_phast/chunks/schMedS3h1*_codons.ss > taffy_phast/schMedS3h1_1st_codons.ss


phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS \
--out-root taffy_phast/schMedS3h1_all_4d.ss taffy_phast/schMedS3h1_all_4d.ss 2> ${OUTNAME}_norm_4d_phylofit.log

# conserved model with the first position
phyloFit --tree "(((schMedS3h1,schMedS3h2),schPol2),(schLug1, schNov1))" --msa-format SS --do-cats 1 \
        --out-root taffy_phast/schMedS3h1_1st_codons taffy_phast/schMedS3h1_1st_codons.ss 2> taffy_phast/schMedS3h1_norm_codons_phylofit.log