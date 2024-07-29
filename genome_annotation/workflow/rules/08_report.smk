rule report_quast:
    input:
        rules.quast.output
    output:
        report(directory("08_report/02_01_quast"), category = C02, subcategory = "01_quast", htmlindex="icarus.html")
    shell:
        """
        mkdir -p 08_report/02_01_quast/
        cp -r 02_general/01_quast/* {output} 
        """

rule report_kraken2:
    input:
        rules.kraken2.output
    output:
        report("08_report/02_02_kraken2/" + ASSEMBLY_ID + ".kraken2.krona.html", 
        caption = "../report/02_02_kraken2.rst", category = C02, subcategory = "02_kraken2")
    conda:
        "../envs/krona.yaml"
    shell:
        """
        ktUpdateTaxonomy.sh
        ktImportTaxonomy -q 2 -t 3 02_general/02_kraken2/kraken2 -o {output}
        """

rule report_blobtools:
    input:
        rules.blobplot.output
    output:
        report([ "08_report/02_04_blobtools/" + ASSEMBLY_ID + i for i in ["_blobplot.stats.txt", "_blobplot.cov0.png", "_blobplot.read_cov.png" ]],
        caption = "../report/02_04_blobtools.rst", category = C02, subcategory = "04_blobtools")
    params:
        indir = "02_general/blobplot/",
        prefix = "08_report/02_04_blobtools/" + ASSEMBLY_ID
    shell:
        """
        cp {params.indir}*.blobDB.json.bestsum.phylum.p8.span.100.blobplot.cov0.png {params.prefix}_blobplot.cov0.png
        cp {params.indir}*.blobDB.json.bestsum.phylum.p8.span.100.blobplot.read_cov.cov0.png {params.prefix}_blobplot.read_cov.png
        cp {params.indir}*.blobDB.json.bestsum.phylum.p8.span.100.blobplot.stats.txt {params.prefix}_blobplot.stats.txt
        """

rule report_busco:
    input:
        rules.busco.output
    output:
        report(["08_report/03_01_busco/" + i for i in ["busco_figure.png", "busco_short_summary.txt"]],
        caption= "../report/03_01_busco.rst", category = C02, subcategory = "01_busco")
    params:
        script = "workflow/scripts/generate_plot.py",
        busco_dir =  "03_completeness/01_busco/busco/",
        out_dir = "08_report/03_01_busco/"
    conda:
        "../envs/busco_plot.yaml"
    shell:
        """
        cp -r {params.busco_dir}*.txt {params.out_dir} &&
        python3 {params.script} -wd {params.out_dir} &&
        cp {params.out_dir}*txt {params.out_dir}busco_short_summary.txt
        """

#rule report_gmap:
#    input:
#        expand("03_completeness/done__03_gmap_{tx}", tx=config['tx_assembly'])


rule ccs_coverage:
    input:
        rules.map_ccs_genome.output


rule r08_report_done:
    input:
       rules.report_quast.output,
       rules.report_kraken2.output,
       #rules.report_blobtools.output,
       rules.report_busco.output
    output:
        touch("08_report/done__08_report")
