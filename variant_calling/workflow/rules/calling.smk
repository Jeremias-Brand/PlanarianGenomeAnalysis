

#TODO reporting structure
onsuccess:
    print("Workflow finished, no error!"),
    shell("mail -s 'snakemake finished OK' jbrand@mpibpc.mpg.de < {log}")

onerror:
    print("An error occurred"),
    shell("mail -s 'snakemake finished NOT OK' jbrand@mpibpc.mpg.de < {log}")

ref_dir = config["ref_dir"]
ref_ext = config["ref_ext"]


vardir = ""

SEdir = vardir + "SE/"
PEdir = vardir + "PE/"
LOGdir = os.path.join(workflow.basedir, "logs/log__variation_")

##############################################################
#   INDEX THE REFERENCE GENOME
############################################################## 

rule index_ref_genome:
    input:
        ref_dir + "{GE}" + ref_ext
    output:
    #TODO set this to temp()
        ref_dir + "{GE}" + ref_ext + ".bwt"
    conda:
        "../envs/06_bwa.yaml"
    log:
        LOGdir + "index_ref_genome_{GE}"
    shell:
        """
        bwa index -a bwtsw {input} &> {log}
        """


rule samdict_ref_genome:
    input:
        ref_dir + "{GE}" + ref_ext
    output:
    #TODO set this to temp()
        ref_dir + "{GE}" + ref_ext + ".dict"
    conda:
        "../envs/06_samtools.yaml"
    log:
        LOGdir + "samdict_ref_genome_{GE}"
    shell:
        """
        samtools dict {input} > {output} 2> {log}
        """


rule seqdict_ref_genome:
    input:
        ref_dir + "{GE}" + ref_ext
    output:
    #TODO set this to temp()
        ref_dir + "{GE}" + ref_ext + ".seqdic.done"
    conda:
        "../envs/gatk.yaml"
    log:
        LOGdir + "seqdict_ref_genome_{GE}" 
    shell:
        """
        gatk CreateSequenceDictionary -R {input} 2> {log} && touch {output}
        """


rule faidx_ref_genome:
    input:
        ref_dir + "{GE}" + ref_ext
    output:
        ref_dir + "{GE}" + ref_ext + ".fai"
    conda:
        "../envs/06_samtools.yaml"
    log:
        LOGdir + "faidx_ref_genome_{GE}"
    shell:
        """
        samtools faidx {input} > {output} 2> {log}
        """




###############################################################################

# QC PE

###############################################################################


# QC

rule PE_fastqc_raw:
    input:
        chck = PEdir + "01_raw/{PE}_validate.done",
        r1 = PEdir + "01_raw/{PE}_1.fastq.gz",
        r2 = PEdir + "01_raw/{PE}_2.fastq.gz"
    output:
        touch(PEdir + "01_raw/{PE}_qc.done")
    params:
        prefix = PEdir + "QC"
    conda:
        "../envs/06_fastqc.yaml"
    log:
        LOGdir + "PE_fastqc_raw_{PE}"
    shell:
        """
        mkdir -p {params.prefix}
        fastqc {input.r1} {input.r2} --outdir {params.prefix} &> {log} 
        """

rule PE_fastqc_trim:
    input:
        r1 = PEdir + "02_trim/{PE}_1.fastq.gz",
        r2 = PEdir + "02_trim/{PE}_2.fastq.gz",
    output:
        touch(PEdir + "02_trim/{PE}_qc.done")
    params:
        prefix = PEdir + "QC"
    conda:
        "../envs/06_fastqc.yaml"
    log:
        LOGdir + "PE_fastqc_trim_{PE}"
    shell:
        """
        mkdir -p {params.prefix}
        fastqc {input.r1} {input.r2} --outdir {params.prefix}
        """

rule PE_bamstats:
    input:
        PEdir + "03_bam/{GE}__{PE}_sorted_markdup_rg.bam"
    output:
        PEdir + "QC/{GE}__{PE}_bamtools.stats"
    conda:
        "../envs/06_bamtools.yaml"
    shell:
        """
        bamtools stats -in {input} | grep -v "*" > {output}
        """

rule PE_deeptools_coverage:
    input:
        expand(PEdir + "03_bam/{{GE}}__{PE}_aligned_duplicates_marked_sorted.bam", PE = config["PE"])
    output:
        plot = PEdir + "QC/deeptools_coverage.pdf",
        raw = PEdir + "QC/deeptools_coverage.raw"
    conda:
        "../envs/06_deeptools.yaml"
    log:
        LOGdir + "_PE_deeptools_coverage"
    threads: 10
    shell:
        """
        plotCoverage -b {input} \
         --plotFile {output.plot} \
         --numberOfSamples 2000000 \
         --outRawCounts {output.raw} \
         --numberOfProcessors {threads} &> {log}
        """


rule PE_fb_deeptools_coverage:
    input:
        expand(PEdir + "03_bam/{{GE}}__{PE}_sorted_markdup_rg.bam", PE = config['PE_fb_selection'])
    output:
        plot = PEdir + "QC/{GE}__PE_fb_deeptools_coverage.pdf",
        raw = PEdir + "QC/{GE}__PE_fb_deeptools_coverage.raw"
    conda:
        "../envs/06_deeptools.yaml"
    log:
        LOGdir + "{GE}__PE_fb_deeptools_coverage"
    threads: 10
    shell:
        """
        plotCoverage -b {input} \
         --plotFile {output.plot} \
         --numberOfSamples 2000000 \
         --outRawCounts {output.raw} \
         --numberOfProcessors {threads} &> {log}
        """


###########################
# Contamination screen
###########################

rule kraken2_index:
    input:
    output:
        touch("logs/done__kraken2_index")
    params:
        "--use-ftp --db tmp/kraken2_index --download-library"
    threads: 1
    log:
        "logs/log__kraken2_index"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        mkdir -p tmp/kraken2_index && cd tmp/kraken2_index &&
        wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz &&
        tar -xvf k2_standard_20201202.tar.gz
        """


rule PE_kraken2:
    input:
        r1 = PEdir + "02_trim/{PE}_1.fastq.gz",
        r2 = PEdir + "02_trim/{PE}_2.fastq.gz",
        k2 = rules.kraken2_index.output
    output:
        touch("logs/done__kraken2_{PE}"),
        f = PEdir + "QC/{PE}.kraken"
    params:
        args = "--db tmp/kraken2_index --confidence 0.5 --minimum-base-quality 20 --use-names --gzip-compressed"
    threads:
        config['threads_kraken2']
    log:
        "logs/log__kraken2_{PE}"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 {params.args} \
        --threads {threads} \
        --report {output.f}.report \
        --paired {input.r1} {input.r2} \
        --output {output.f}  > {log} 2>&1
        """



rule PE_bcftools_stats_single:
    input:
        PEdir + "06_vcf/single_all_raw.vcf.gz"
    output:
        PEdir + "QC/bcftools_stats_single_all_raw_vcf.txt"
    conda:
        "../envs/06_bcftools.yaml"
    log:
        LOGdir + "PE_bcftools_stats_single"
    shell:
        """
        bcftools stats {input} > {output} &> {log}
        """

rule PE_bcftools_stats_joint:
    input:
        PEdir + "06_vcf/joint_all_raw.vcf.gz"
    output:
        PEdir + "QC/bcftools_stats_joint_all_raw_vcf.txt"
    conda:
        "../envs/06_bcftools.yaml"
    log:
        LOGdir + "PE_bcftools_stats_joint"
    shell:
        """
        bcftools stats {input} > {output} &> {log}
        """


rule PE_kat:
    input:
        r1 = PEdir + "02_trim/{PE}_1.fastq.gz",
        r2 = PEdir + "02_trim/{PE}_2.fastq.gz"
    output:
        PEdir + "06_kmer/{PE}_kat.hist"
        # jf has will be called:
        # PEdir + "06_kmer/{PE}_kat.hist-hash.jf21"
    params:
        "-m 21 -H 1000000000 -d "
    conda:
        # some issue with the conda recipie had to use 2.4.1
        "../envs/kat.yaml"
    log:
        LOGdir + "PE_kat_{PE}"
    threads: config['threads_kat']
    shell:
        """
        kat hist {params} -t {threads} -o {output} {input.r1} {input.r2} > {log}
        tail -n +7 {output} > {output}_jf.hist
        """


rule PE_kat_gcp:
    input:
        r1 = PEdir + "02_trim/{PE}_1.fastq.gz",
        r2 = PEdir + "02_trim/{PE}_2.fastq.gz"
    output:
        touch(PEdir + "QC/done__PE_kat_gcp_{PE}")
    params:
        prefix = PEdir + "QC/kat_gcp_",
        par = "-m 27 -H 1000000000"
    conda:
        "../envs/kat.yaml"
    log:
        LOGdir + "PE_kat_gcp_{PE}"
    threads: config['threads_kat']
    shell:
        """
        kat gcp -o {params.prefix}{wildcards.PE} {params.par} -t {threads} {input.r1} {input.r2} > {log}
        """


rule PE_genomescope:
    input:
        rules.PE_kat.output
    output:
        touch(directory(PEdir + "07_genomescope/done__genomescope_{PE}")),
    params:
        script = "workflow/scripts/genomescope.R",
        read_length = config["PE_read_length"]
    log:
        LOGdir + "genomescope_{PE}"
    threads: config['threads_genomescope']
    conda:
        "../envs/r.yaml"
    shell:
        """
        outdir={output}
        prefix={wildcards.PE}_
        Rscript {params.script} {input}_jf.hist 21 {params.read_length} $outdir 1000 1 > {log} 2>&1
        mv $outdir/summary.txt $outdir/${{prefix}}summary.txt
        mv $outdir/plot.png $outdir/${{prefix}}plot.png
        mv $outdir/plot.log.png $outdir/${{prefix}}plot.log.png
        # when coverage is low genomescope often does not converge
        # so we only move this file if it was generated
        if [ -f $outdir/model.txt ]; then
            mv $outdir/model.txt $outdir/${{prefix}}model.txt
        fi
        """


rule PE_multiqc:
    input:
        expand(PEdir + "01_raw/{PE}_qc.done", PE=config["PE_fb_selection"]),
        expand(PEdir + "02_trim/{PE}_qc.done", PE=config["PE_fb_selection"]),
        expand(PEdir + "QC/{GE}__{PE}_bamtools.stats", PE=config["PE_fb_selection"], GE=config["PE_fb_ref"]),
        expand(PEdir + "QC/{GE}__PE_fb_deeptools_coverage.pdf", GE=config["PE_fb_ref"]),
        expand(PEdir + "QC/{PE}.kraken", PE=config["PE_fb_selection"]),
        expand(PEdir + "QC/done__PE_kat_gcp_{PE}", PE=config['PE_fb_selection']),
        expand(PEdir + "07_genomescope/done__genomescope_{PE}", PE=config['PE_fb_selection'])
    output:
        PEdir + "PE_multiqc.html"
    params:
        d = PEdir,
        name = "PE_multiqc.html"
    conda:
        "../envs/06_multiqc.yaml"
    log:
        LOGdir + "PE_multiqc"
    shell:
        """
        cd {params.d}
        multiqc --force --filename {params.name} QC/
        """


###############################################################################

# freebayes align PE

###############################################################################

rule PE_prefetch:
    output:
        touch(PEdir + "01_raw/.prefetch/sra/{PE}.sra")
    params:
        args = "{PE} --max-size 50GB",
        dir = PEdir + "01_raw/"
    log:
        LOGdir + "PE_prefetch_{PE}"
    conda:
        "../envs/sra-tools.yaml"
    shell:
        """
        mkdir -p {params.dir}
        cd {params.dir}
        prefetch {params.args} > {log} 2>&1
        """


rule PE_fastqdump:
# TODO check the gzip works
    input:
        PEdir + "01_raw/.prefetch/sra/{PE}.sra"
    output:
        r1 = PEdir + "01_raw/{PE}_1.fastq.gz",
        r2 = PEdir + "01_raw/{PE}_2.fastq.gz"
    params:
        args = "-S -O 01_raw",
        id_srr = "{PE}",
        dir = PEdir + "01_raw/"
    log:
        LOGdir + "PE_fastqdump_{PE}"
    conda:
        "../envs/sra-tools.yaml"
    shell:
        """
        cd {params.dir}
        fasterq-dump {params.args} {params.id_srr} > {log} 2>&1
        gzip *.fastq
        """



rule PE_validate:
    input:
        r1 = PEdir + "01_raw/{PE}_1.fastq.gz",
        r2 = PEdir + "01_raw/{PE}_2.fastq.gz"
    output:
        touch(PEdir + "01_raw/{PE}_validate.done")
    conda:
        "../envs/06_validate.yaml"
    params:
        args = "--log_level info"
    log:
        LOGdir + "PE_validate_{PE}"
    shell:
        """
        biopet-validatefastq {params.args} --fastq1 {input.r1} --fastq2 {input.r2} &> {log}
        """


rule PE_trimmomatic:
    input:
        PEdir + "01_raw/{PE}_validate.done",
        r1 = PEdir + "01_raw/{PE}_1.fastq.gz",
        r2 = PEdir + "01_raw/{PE}_2.fastq.gz"
    output:
        r1 = PEdir + "02_trim/{PE}_1.fastq.gz",
        r2 = PEdir + "02_trim/{PE}_2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unp = PEdir + "02_trim/{PE}_1_unpaired.fastq.gz",
        r2_unp = PEdir + "02_trim/{PE}_2_unpaired.fastq.gz"
    conda:
        "../envs/06_trimmomatic.yaml"
    log:
        LOGdir + "PE_trimmomatic_{PE}"
    params:
        # list of trimmers (see manual)
        trimmer= config["trimmomatic_args"],
    threads:
        config["threads_trimmomatic"]
    shell:
        """
        trimmomatic PE -threads {threads} \
        {input.r1} {input.r2} \
        {output.r1} {output.r1_unp} \
        {output.r2} {output.r2_unp} \
        {params.trimmer} &> {log}
        """

rule PE_fb_bwa:
    input:
        r1 = PEdir + "02_trim/{PE}_1.fastq.gz",
        r2 = PEdir + "02_trim/{PE}_2.fastq.gz",
        ref_genome = ref_dir + "{GE}" + ref_ext,
        ref_genome_index = ref_dir + "{GE}" + ref_ext + ".bwt",
        fai = ref_dir + "{GE}" + ref_ext + ".fai"
    output:
        temp(PEdir + "03_bam/{GE}__{PE}_aligned_rg.bam")
    conda:
        "../envs/freebayes.yaml"
    threads: 30
    log:
        LOGdir + "PE_fb_bwa_{GE}__{PE}"
    shell:
        """
        bwa mem -t {threads} \
        {input.ref_genome} \
        {input.r1} \
        {input.r2} \
        -K 100000000 \
        -R  '@RG\\tID:{wildcards.PE}\\tLB:{wildcards.PE}\\tSM:{wildcards.PE}\\tPL:ILLUMINA' \
        -v 0 2> {log} | \
        samtools view -1 - > {output} && [[ -s {output} ]]
        """


rule PE_fb_sort:
    input:
        rules.PE_fb_bwa.output
    output:
        temp(PEdir + "03_bam/{GE}__{PE}_sorted_rg.bam")
    log:
        LOGdir + "PE_fb_sort_{GE}__{PE}"
    threads: 8
    wrapper:
        "v1.3.2/bio/sambamba/sort"


rule PE_fb_markdup:
    input:
        rules.PE_fb_sort.output
    output:
        PEdir + "03_bam/{GE}__{PE}_sorted_markdup_rg.bam"
    params:
        extra="-r"  # optional parameters
    log:
        LOGdir + "PE_fb_markdup_{GE}__{PE}"
    threads: 8
    wrapper:
        "v1.3.2/bio/sambamba/markdup"


rule PE_freebayes_gvcf:
    input:
        ref = ref_dir + "{GE}" + ref_ext,
        ref_genome_index = ref_dir + "{GE}" + ref_ext + ".bwt",
        # we want a combination of wildcard and expand. we do this using {{}}
        samples=expand(PEdir + "03_bam/{{GE}}__{PE}_sorted_markdup_rg.bam", PE = config['PE_fb_selection'])
    output:
        PEdir + "04_gvcf/{GE}__" + "_".join(config['PE_fb_selection']) +  "_freebayes.g.bcf"  # either .vcf or .bcf
    log:
        LOGdir + "PE_freebayes_gvcf_{GE}__" + "_".join(config['PE_fb_selection'])   # either .vcf or .bcf
    params:
        extra="--gvcf -g 1000",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: config['threads_freebayes']
    conda: 
        "../envs/freebayes.yaml"
    wrapper:
        "v1.2.0/bio/freebayes"


rule PE_freebayes_gvcf_all:
    input: 
        expand(PEdir + "04_gvcf/{GE}__" + "_".join(config['PE_fb_selection']) +  "_freebayes.g.vcf", GE = config['PE_fb_ref'])  # either .vcf or .bcf
    output:
        touch("done__freebayes_gvcf_all")


