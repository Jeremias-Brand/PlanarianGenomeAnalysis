# we run orthofinder via genespace and include BraLan, NemVec and the parasites
# activate env
conda activate genespace
# run GENESPACE
Rscript --vanilla ./scripts/orthofinder_via_genespace.R
# copy result file
cp ./genespace_parasites_metazoa/results/combBed.txt orthofinder/combBed.txt 
# plot and perform statistical tests
Rscript --vanilla ./scripts/orthofinder_significance_tests_synteny_og.R
Rscript --vanilla ./scripts/orthofinder_dotplots_from_og.R