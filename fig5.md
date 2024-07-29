
# ODP

- Run `scripts/setup_odp.sh`.
- Prepare reference sequences in the folder structure expected by ODP. 
- Prepare additional files and references using `scripts/ODP.sh`

Run ODP using the config `scripts/ODP_config.yaml`

```bash
snakemake -r -p --snakefile odp/scripts/odp -n
```
## Replotting to make Fig5 and supporting figures

We replot the figures using `scripts/fig5_odp_replot.R`.
