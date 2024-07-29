# Figure 4a, 4c

EDTA and other annotation pipelines are in the `genome_annotation` directory. 
Once the annotation is done the plots can be produced using `scripts/fig4_EDTA.R`.

# Figure 4b

Run `scripts/setup_genespace.sh`, `scripts/setup_genespace_directory.sh`, `scripts/genespace.R`, and `scripts/genespace_plotting.R`.

# Figure 4d

Run orthofinder via genespace `scripts/species_tree_orthofinder.R`.

Once done we realign the single copy genes.
Then we concatenate and infer the phylogeny.

```bash
cd ./species_tree/phylogeny

mamba create -n phy mafft iqtree

# mafft E-INS-i
mkdir -p aln
for a in $( ls Single_Copy_Orthologue_Sequences/*.fa );do
name=$( basename $a .fa )
echo "mafft --ep  0 --genafpair --maxiterate 1000 $a > aln/${name}.aln"
done >> mafft.cmd
parallel -a mafft.cmd -j 144 --joblog mafft.joblog

mkdir -p  renamed_aln 
# making uniform names between files
for a in $( ls aln/*.aln );do
name=$( basename $a .aln )
sed -E '/>/s/_?([0-9]+-)?[0-9.T]+$//' $a | sed -E '/>/s/>Tm.G/>TmG/'  > renamed_aln/${name}.aln &
done 

iqtree -T 144 --safe -m MFP -p renamed_aln/ -B 1000 -bnni -alrt 1000 --prefix species_tree
```

# Figure 4e

Execute the snakemake workflow in `busco_synteny`.
For the supporting analysis with orthofinder-based orthologs first run the analysis for Figure 4b, then run `scripts/orthofinder_dotplots_from_og.R` and `scripts/orthofinder_significance_tests_synteny_og.R`.
