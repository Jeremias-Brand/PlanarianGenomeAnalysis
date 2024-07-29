rule busco:
    input:
        fas = rules.rename_genome_assembly.output
    output:
        touch("03_completeness/01_busco/done__busco")
    log:
        os.path.join(workflow.basedir, "logs/log__03_completeness_01_busco")
    params:
        args = "-m genome -o busco -l eukaryota",
        outdir = "03_completeness/01_busco"
    threads:
        config['threads_busco']
    conda:
        "../envs/busco.yaml"
    shell:
        """
        mkdir -p {params.outdir} && 
        cd {params.outdir} && 
        busco -i ../../{input.fas} -c {threads} {params.args} > {log} 2>&1
        """

rule gmap_index:
    input:
        config['genome_fasta_file']
    output:
        touch("logs/done__03_gmap_index")
    log:
        "logs/log__03_gmap_index"
    conda:
        "../envs/gmap.yaml"
    shell:
        """
        gmap_build -D tmp/ -d gmap_index {input} > {log} 2>&1
        """


rule gmap:
    input:
        fasta = "lib/tx/{tx}.PCFL.fasta",
        ix = rules.gmap_index.output
    output:
        touch("03_completeness/done__03_gmap_{tx}")
    params:
        idx = "-D tmp -d gmap_index",
        args = "--chimera-margin=200 --min-trimmed-coverage=0.6 --min-identity=0.6 -F -Y -f gff3_gene",
        out = "03_completeness/02_gmap/{tx}"
    threads: config['threads_gmap']
    log:
        "logs/log__03_gmap_{tx}"
    conda:
        "../envs/gmap.yaml"
    shell:
        """
        mkdir -p 03_completeness/02_gmap
        gmap {params.idx} {params.args} --split-output={params.out} -t {threads} {input.fasta} > {log} 2>&1
        """

rule r03_completeness_done:
    input:
        rules.busco.output,
        expand("03_completeness/done__03_gmap_{tx}", tx=config['tx_assembly'])
    output:
        touch("03_completeness/done__03_completeness")
