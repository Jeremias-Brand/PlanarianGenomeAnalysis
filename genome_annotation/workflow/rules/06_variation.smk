rule minimap2:
    input:
        genome = rules.rename_genome_assembly.output,
        ccs = rules.merge_ccs.output
    output:
        "06_variation/02_sv/output.bam"
    log:
        "logs/log__06_variation_02_sv_mm2"
    conda:
        "../envs/mm2.yaml"
    threads:
        config['threads_mm2']
    shell:
        """
        minimap2 -t {threads} -ax asm20 --MD -R @RG'\\'tID:{config[mm2_id]}'\\'tSM:{config[mm2_sample]} {input.genome} {input.ccs} | samtools view -bh - |samtools sort -@4 - > {output}
        """

rule sniffles:
    input:
        rules.minimap2.output
    output:
        "06_variation/02_sv/sniffles.vcf"
    log:
        "logs/log__06_variation_02_sv_sniffles"
    conda:
        "../envs/sniffles.yaml"
    threads:
        config['threads_sniffles']
    shell:
        """
        sniffles -s 5 -m {input} --report_BND -v {output} -t {threads} > {log} 2>&1
        """

rule svim:
    input:
        fasta = rules.rename_genome_assembly.output,
        bam = rules.minimap2.output
    output:
        touch("06_variation/02_sv/done__svim")
    params:
        "06_variation/02_sv/svim"
    log:
        "logs/log__06_variation_02_sv_svim"
    conda:
        "../envs/svim.yaml"
    threads:
        config['threads_svim']
    shell:
        """
        svim alignment {params} {input.bam} {input.fasta} > {log} 2>&1
        """


rule r06_variation_done:
    input:
        rules.minimap2.output,
        rules.sniffles.output,
        rules.svim.output
    output:
        touch("06_variation/done__06_variation")

