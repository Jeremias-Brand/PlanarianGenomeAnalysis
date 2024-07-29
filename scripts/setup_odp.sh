# setup mamba env
mamba create -y -n odp snakemake matplotlib pandas numpy seaborn diamond blast hmmer networkx mscorefonts
conda activate odp

# download the repo and make
git clone https://github.com/conchoecia/odp.git
cd odp 

## there is an error in the make file because it does not specify cores
# adjust depending on available cores:
sed -i 's/snakemake/snakemake --cores 144/' Makefile
make

exit
## To address this issue https://github.com/conchoecia/odp/issues/34
## We make a change to the odp main file ./odp/scripts/odp

## Fix from github issue 34
'''
from matplotlib import font_manager


font_dirs = ["/home/jbrand/anaconda3/envs/odp/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/pdfcorefonts"]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
matplotlib.use('SVG')
'''

## We also change one rule:

rule filtered_D_FET_rbh:
    """
    This rule filters the RBH file with D and FET calculated to generate an RBH with just the FET significant hits.
    We also should assign a group name to these significant groups.
    """
    input:
        rbh = config["tool"] + "/step1-rbh/{analysis}_reciprocal_best_hits.D.FET.rbh"
    output:
        filt = config["tool"] + "/step1-rbh-filtered/{analysis}_reciprocal_best_hits.D.FET.filt.rbh"
    threads: 1
    run:
        # Read in the rbh file
        df = pd.read_csv(input.rbh, delimiter="\t")
        sp_order = wildcards.analysis.split("_")

        # Filter out the rows where the whole_FET is less than or equal to 0.05
        df = df[df["whole_FET"] <= 0.05]

        # Groupby the x and y scaf column values
        scaf_cols = ["{}_scaf".format(x) for x in sp_order]
        df["color"] = "#000000"

        gb = df.groupby(scaf_cols)

        # Iterate through the groups and assign a group name
        for name, group in gb:
            gene_group = "{}_{}".format(name[0], name[1])
            # Change the gene_group value of this group in the group, not the df
            group.loc[:, "gene_group"] = gene_group
            group.loc[:, "color"] = grc()

        # Unwrap the groups back into a df
        df = gb.first().reset_index(drop=True)

        # Filter the df based on whole_FET values
        df = df[df["whole_FET"] <= 0.05]

        # Save the filtered df to a CSV file
        df.to_csv(output.filt, sep="\t", header=True, index=None)

