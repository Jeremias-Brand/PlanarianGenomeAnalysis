rule fetch_sra_all:
    input:
        #expand("01_raw/.prefetch/sra/{srr}.sra", srr=config['srr'])
        expand("lib/sra/done__fastqdump_{srr}", srr=config['fetch_srr'])
    output:
        touch("lib/sra/done__fetch_sra_all")

rule prefetch:
    output:
        touch("lib/sra/done__{srr}.sra")
    params:
        "{srr} --max-size 50GB -O lib/sra/"
    log:
        "lib/sra/log__prefetch_{srr}"
    conda:
        "../envs/sra-tools.yaml"
    shell:
        """
        prefetch {params} > {log} 2>&1
        """

rule fastqdump:
    input:
        "lib/sra/done__{srr}.sra"
    output:
        touch("lib/sra/done__fastqdump_{srr}")
    params:
        args = "-S -O lib/sra/ -t lib/sra",
        id_srr = "{srr}"
    threads:
        10
    log:
        "lib/sra/log__fastqdump_{srr}"
    conda:
        "../envs/sra-tools.yaml"
    shell:
        """
        fasterq-dump --threads {threads} {params.args} {params.id_srr} > {log} 2>&1
        """

