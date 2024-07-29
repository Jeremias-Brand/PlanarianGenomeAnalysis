# 2022-gqc

Slim genome QC project.
Builds on the previous genome_check pipeline.

## Setup
The workflow uses [conda](https://www.anaconda.com/) and [snakemake](https://snakemake.readthedocs.io/en/stable/).
To get started download and install conda.
We also recommend that you use [mamba](https://github.com/mamba-org/mamba) a faster way to install conda packages.

`conda install mamba -n base -c conda-forge`

Then create a environment for snakemake and activate it:
```
mamba create -y -n snakemake snakemake mamba
conda activate snakemake
```

Edit the config/config.yaml file to tell the workflow the assemblies you want to check.
For each assembly the tools specified below will be run.

Once you are set up you can run the complete workflow like this:
```
snakemake --use-conda --cores 20
```

Or if you prefere to only run parts of the workflow, then you can edit the Snakemake file and comment out certain lines in the all rule. 
