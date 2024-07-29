rule busco:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_busco/done__busco_{GE}")
    log:
        os.path.join(workflow.basedir, "01_busco/log__busco_{GE}")
    params:
        args = "-m genome -o busco -l metazoa",
        outdir = "01_busco"
    threads:
        config['threads_busco']
    conda:
        "../envs/busco.yaml"
    shell:
        """
        d={params.outdir}/{wildcards.GE}
        mkdir -p $d && 
        cd $d && 
        busco -f -i ../../{input.fas} -c {threads} {params.args} > {log} 2>&1
        """

rule all_busco:
    input:
        expand("01_busco/done__busco_{GE}", GE=config['genomes'])
    output:
        touch("done__busco")
    shell:
        """
        find 01_busco -name '*down*' -exec rm -rdf {} +
        """
