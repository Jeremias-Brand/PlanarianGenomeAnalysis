rule md5sum_orig:
    input:
        config['genome_fasta_file']
    output:
        "01_naming/orig.md5sum"
    threads: 1
    shell:
        """
        md5sum {input} > {output}
        """

rule rename_genome_assembly:
    input:
        fasta = config['genome_fasta_file'],
        md5sum = rules.md5sum_orig.output
    output:
        "01_naming/"+config['ucsc_id']+a"
    params:
        ucsc_id = config['ucsc_id'],
        name = config['name']
    log:
        "logs/log__01_naming__rename_genome_assembly"
    threads: 1
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python workflow/scripts/rename_fasta.py {input.fasta} {output} {params.ucsc_id} {params.name} > {log} 2>&1
        """

rule md5sum_renamed:
    input:
        rules.rename_genome_assembly.output
    output:
        report("01_naming/"+config['ucsc_id']+".fasta.md5sum",caption="../report/01_naming.rst",category="01 Naming",subcategory="03 md5sum renamed")
    threads: 1
    shell:
        """
        md5sum {input} > {output}
        """

rule r01_naming_done:
    input:
        rules.md5sum_orig.output,
        rules.md5sum_renamed.output
    output:
        touch("01_naming/done__01_naming")
