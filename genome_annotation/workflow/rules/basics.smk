rule gunzip:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT'] + ".gz"
    output:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    shell:
        """
        gunzip {input}
        """


rule gzip:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch(config['genomeDIR'] + "done__gzip_{GE}")
    shell:
        """
        gzip {input}
        """

rule all_gzip:
    input:
        expand(config['genomeDIR'] + "{GE}" + config['genomeEXT'], GE=config['genomes'])

rule assembly_stats:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        "01_assembly_stats/{GE}.stats"
    params:
        "-t"
    conda:
        "../envs/assembly_stats.yaml"
    shell:
        """
        assembly-stats {params} {input} > {output}
        """


rule quast:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        directory("01_quast/{GE}")
    log:
        "01_quast/log__quast_{GE}"
    params:
        "01_quast"
    threads:
        config['threads_quast']
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast -o {output} -t {threads} {input} > {log} 2>&1
        """


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
        busco -i ../../{input.fas} -c {threads} {params.args} > {log} 2>&1
        """


rule EDTA:
    input:
        fasta = os.path.join(workflow.basedir, config['EDTADIR'] + '{GE}' + config['genomeEXT']),
        cds = os.path.join(workflow.basedir, config['EDTADIR'] + '{GE}.cds') 
    output:
        touch("01_EDTA/done__{GE}")
    params:
        dir = '01_EDTA',
	# remove the force with large targets
        # arg = "--species others --step all --sensitive 1 -anno 1",
        arg = "--species others --step all --sensitive 0 -anno 1",
        curatedlib = os.path.join(workflow.basedir, config['curatedlib'])
    log:
        os.path.join(workflow.basedir, 'log/EDTA__{GE}')
    conda:
        "../envs/EDTA.yaml"
    threads:
        config['threads_EDTA']
    shell:
        """
        cd {params.dir}
        mkdir -p {wildcards.GE}
        cd {wildcards.GE}
        EDTA.pl --genome {input.fasta} -t {threads} \
        --cds {input.cds} \
        --curatedlib {params.curatedlib} \
        {params.arg} &> {log}
        """


rule rm:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_rm/done__rm_{GE}")
    params:
        outdir = "01_rm",
        lib = os.path.join(workflow.basedir, config['repeatmasker_lib']),
        args = "-a -xsmall -no_is -s -e rmblast -dir ./ -gff"
    log:
        "logs/01_rm_{GE}"
    threads: 
        config['threads_repeatmasker']
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        mkdir -p {params.outdir} && cd {params.outdir} && 
        RepeatMasker -pa {threads} {params.args} -lib {params.lib} ../{input} > ../{log} 2>&1
        """

# Mappability with gem
rule gem_index:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        idx = "01_mappability/{GE}_gem_index.gem"
    params:
        outdir = "01_mappability",
        idx = "01_mappability/{GE}_gem_index",
        gem_bin = os.path.join(workflow.basedir, "workflow/scripts/bin")
    log:
        "logs/log__01_mappability_{GE}_gem_index"
    threads:
        config['threads_gem_index']
    shell:
        """
        export PATH={params.gem_bin}:$PATH
        mkdir -p {params.outdir} && 
        workflow/scripts/bin/gem-indexer -T {threads} -c dna -i {input} -o {params.idx} > {log} 2>&1
        """

rule gem:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        idx = rules.gem_index.output
    output:
        touch("01_mappability/done__{GE}_gem")
    params:
        idx = "01_mappability/{GE}_gem_index",
        out = "01_mappability/{GE}_gem_out",
        gem_bin = os.path.join(workflow.basedir, "workflow/scripts/bin")
    log:
        "logs/log__01_mappability_{GE}_gem"
    threads:
        config['threads_gem']
    shell:
        """
        export PATH={params.gem_bin}:$PATH
        gem-mappability -T {threads} -I {input.idx} -l 100 -o {params.out} > {log} 2>&1
        """

 
rule gem_to_bw:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        gem = rules.gem.output
    output:
        touch("01_mappability/done__{GE}_bw")
    params:
        idx = "01_mappability/{GE}_gem_index",
        gem_out = rules.gem.params.out,
        gem_bin = os.path.join(workflow.basedir, "workflow/scripts/bin")
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/log__01_mappability_{GE}_bw"
    shell:
        """
        export PATH={params.gem_bin}:$PATH
        samtools faidx {input.genome} 
        cat {input.genome}.fai | awk '{{print$1"\t"$2}}' > {input.genome}.chrsizes
        gem-2-wig -I {params.idx}.gem -i {params.gem_out}.mappability -o {params.gem_out} 2>> {log} 
        cat {params.gem_out}.wig | awk -v OFS="\t" '{{print$1,$2,$5}}' > {params.gem_out}1.wig
        wigToBigWig {params.gem_out}1.wig {input.genome}.chrsizes {params.gem_out}.bw > {log} 2>&1
        """

###############################################################################

# KMER BASED MAPPABILITY ESTIMATION

###############################################################################rule gem_to_bw:

rule khmer_genome_size:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        expand("01_mappability/{{GE}}_khmer_genome_size_r{N}", N = list(range(25,1025,25)))
    params:
        "01_mappability/{GE}_khmer_genome_size_r"
    conda:
        "../envs/khmer.yaml"
    log:
        "logs/log__{GE}_khmer_genome_size"
    shell:
        """
        for i in `seq 25 25 1000`; do
            unique-kmers.py -k ${{i}} --report {params}${{i}} {input} &> {log}
        done
        """


rule khmer_gather:
    input:
        rules.khmer_genome_size.output
    output:
        "01_mappability/{GE}_khmer_genome_size.report",
    conda:
        "../envs/khmer.yaml"
    log:
        "logs/log__{GE}_khmer_gather"
    shell:
        """
        echo {wildcards.GE}
        for k in {input};do
            echo $k
            head -n2 $k | tail -n1 | awk -v var={wildcards.GE} '{{printf "%s\\t%s\\t%s\\n", var, $2, $1}}' >> {output}
        done
        """


rule khmer_plot:
    input:
        rules.khmer_gather.output
    output:
        "01_mappability/{GE}_khmer_genome_size.pdf"
    params:
        os.path.join(workflow.basedir, "workflow/scripts/plot_khmer_report.R")
    conda:
        "../envs/r.yaml"
    log:
        "logs/log__{GE}_khmer_plot"
    shell:
        """
        Rscript --vanilla {params} --infile {input} --outfile {output} &> {log}
        """


###############################################################################

# BLAST SMALL LIST OF COMPLETE FLATWORM MITOCHONDRIAL GENOMES AGAINST ASSEMBLY

###############################################################################

rule make_genome_blastdb:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_blastn_mito/done__{GE}_make_genome_blastdb")
    params:
        db_name = "{GE}",
        assembly_name = "{GE}" + config['genomeEXT'],
        db_prefix = "{GE}"
    log:
        os.path.join(workflow.basedir, "log/log__{GE}_make_genome_blastdb")
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        mkdir -p tmp/blastDB_genome/
        cp {input} tmp/blastDB_genome/
        cd tmp/blastDB_genome/
        makeblastdb \
        -dbtype nucl \
        -input_type fasta \
        -out {params.db_prefix} \
        -title {params.db_name} \
        -in {params.assembly_name} > {log} 2>&1
        """


rule unpack_mito:
    input:
        mito_db = "workflow/resources/plan_mitogenomes_kraken_db.fa.gz"
    output:
        temp("tmp/mito_query.fa")
    shell:
        """
        zcat {input} > {output}
        """


rule blastn_mito_to_genome:
    input:
        rules.make_genome_blastdb.output,
        assembly = config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        mito_db = rules.unpack_mito.output
    output:
        touch("01_blastn_mito/done__{GE}_blastn_mito"),
        tbl = report("01_blastn_mito/mito_to_{GE}_blast.out")
        #caption = "../report/04_06_blast_mito_to_genome.rst", category = "04_annotation", subcategory = "06_blastn_mito")
    params:
        db_name = "tmp/blastDB_genome/{GE}",
    threads:
        config['threads_blastn']
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        blastn -db {params.db_name} \
        -num_threads {threads} \
        -query {input.mito_db} \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qcovs sstart send evalue bitscore' > {output.tbl}
        echo "contigs in "{wildcards.GE} "  with a 95% query coverage of a mitogenome:" > {output.tbl}.summary
        awk '$9> 95 {{print $2}}' {output.tbl} | uniq >> {output.tbl}.summary
        echo "contigs in "{wildcards.GE} "  with a 99% query coverage of a mitogenome:" >> {output.tbl}.summary
        awk '$9> 99 {{print $2}}' {output.tbl} | uniq >> {output.tbl}.summary
        """


rule satDNA:
    # This is simply a repeatmasker run with a sattelite repeat library
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_satDNA/done__satDNA_{GE}")
    params:
        outdir = "01_satDNA",
        lib = os.path.join(workflow.basedir, config['satDNA_lib']),
        args = "-xsmall -no_is -e ncbi -nolow -dir ./ -gff"
    log:
        "logs/01_satDNA_{GE}"
    threads: 
        config['threads_repeatmasker']
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        mkdir -p {params.outdir} && cd {params.outdir} && 
        RepeatMasker -pa {threads} {params.args} -lib {params.lib} ../{input} > ../{log} 2>&1
        """

rule satDNA_cleanup:
    # https://github.com/kavonrtep/repeat_annotation_pipeline/blob/main/clean_rm_output.R
    input:
        "01_satDNA/done__satDNA_{GE}"
    output:
        "01_satDNA/{GE}_satDNA.gff3"
    params:
        outdir = "01_satDNA",
        script = os.path.join(workflow.basedir, 'workflow/scripts/clean_rm_output.R')
    log:
        "logs/01_satDNA_cleanup_{GE}"
    threads: 
       8 # It's hardcoded... 
    conda:
        "../envs/r-cleanup.yaml"
    shell:
        """
        Rscript --vanilla {params.script} {params.outdir}/{wildcards.GE}.fa.out {output}
        """


rule satDNA_filter:
    input:
       rules.satDNA_cleanup.output
    output:
        o1 = "01_satDNA/{GE}_satDNA_min1k.gff3", 
        o2 = "01_satDNA/{GE}_satDNA_min1k_remaining.gff3"
    params:
        "--test ">=" --size 1000"
    log:
        "logs/01_satDNA_filter_{GE}"
    conda:
        "../envs/agat.yaml"
    shell:
        """
        agat_sp_filter_gene_by_length.pl --gff {input} {params} --out {output.o1}
        """


rule barrnap:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        fa = "01_barrnap/{GE}_rRNA.fa",
        gff = "01_barrnap/{GE}_rRNA.gff",
        fa2 = "01_barrnap/{GE}_rRNA_bac.fa",
        gff2 = "01_barrnap/{GE}_rRNA_bac.gff"
    params:
        kingdom = config['barrnap_kingdom']
    log:
        "logs/01_barrnap_{GE}"
    threads: 
        config['threads_barrnap']
    conda:
        "../envs/barrnap.yaml"
    shell:
        """
        barrnap --kingdom 'euk' --threads {threads} \
         --outseq {output.fa} < {input} > {output.gff} 2> {log}
        barrnap --kingdom 'bac' --threads {threads} \
         --outseq {output.fa2} < {input} > {output.gff2} 2> {log}
        """


###############################################################################

## SRF

###############################################################################

rule setup_srf:
    output:
        touch('01_srf/done__setup_srf')
    params:
        repo = "https://github.com/lh3/srf"
    log:
        "logs/01_setup_srf"
    threads: 1
    conda:
        "../envs/srf.yaml"
    shell:
        """
        cd 01_srf
        git clone {params.repo}
        cd srf && make && chmod +777 ./srf
        """


rule srf_kmc:
    input:
        rules.setup_srf.output, 
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        t = "01_srf/{GE}.txt"
    params:
        tmp = "./tmp",
        kmc = '-fm -k151 -ci100 -cs100000'
    log:
        "logs/01_srf_kmc_{GE}"
    threads: 
        config['threads_kmc']
    conda:
        "../envs/srf.yaml"
    shell:
        """
        mkdir -p 01_srf/tmp_{wildcards.GE}
        kmc  {params.kmc} -t{threads}  {input.genome} 01_srf/{wildcards.GE}.kmc 01_srf/tmp_{wildcards.GE} 2> {log}
        kmc_dump 01_srf/{wildcards.GE}.kmc {output.t} 2>> {log}
        """


rule srf_srf:
    input:
        rules.srf_kmc.output.t
    output:
        "01_srf/{GE}.srf.fa"
    params:
        script = "01_srf/srf/srf"
    log:
        "logs/01_srf_srf_{GE}"
    threads: 1
    conda:
        "../envs/srf.yaml"
    shell:
        """
        {params.script} -p {wildcards.GE} {input} > {output} 2> {log}
        """


rule srf_mm2:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        srf = rules.srf_srf.output
    output:
        paf = "01_srf/{GE}_srf.paf"
    params:
        script = "01_srf/srf/srf",
        util = "01_srf/srf/srfutils.js",
        mm2 = '-c -N1000000 -f1000 -r100,100'
    log:
        "logs/01_srf_mm2_{GE}"
    threads: 
        config['threads_minimap2']
    conda:
        "../envs/srf.yaml"
    shell:
        """
        minimap2 -t {threads} {params.mm2} <( {params.util} enlong {input.srf} ) {input.genome} > {output.paf} 2> {log}
        """


rule srf_util:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        paf = rules.srf_mm2.output
    output:
        bed = "01_srf/{GE}_srf.bed",
        len = "01_srf/{GE}_srf.len",
        mask = "01_srf/{GE}_srf_masked.fa.gz"
    params:
        script = "01_srf/srf/srf",
        util = "01_srf/srf/srfutils.js",
        mm2 = '-c -N1000000 -f1000 -r100,100'
    log:
        "logs/01_srf_util_{GE}"
    threads: 
        config['threads_minimap2']
    conda:
        "../envs/srf.yaml"
    shell:
        """
        # translate alignment into bed file
        {params.util} paf2bed {input.paf} > {output.bed} 2>> {log} 

        # get stats on the bed file
        {params.util} bed2abun {output.bed} > {output.len} 2>> {log} 

        bedtools maskfasta -fi {input.genome} -bed {output.bed} -fo 01_srf/{wildcards.GE}_srf_masked.fa -soft && gzip 01_srf/{wildcards.GE}_srf_masked.fa
        """



###############################################################################

# GET K-MER SPECTRUM AND RUN GENOMESCOPE

###############################################################################

rule jellyfish_ccs:
    input:
        fof = "config/{GE}.FoF",  # Assuming each GE has a corresponding FoF
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("02_general/{GE}/done__03_jellyfish"),
        jf_k21 = "02_general/{GE}/03_genomescope/genomic_ccs_k21.jf",
        hist = "02_general/{GE}/03_genomescope/{GE}.hist"
    params:
        # No need for the reads parameter here since we'll be reading from the FoF
    threads: 
        config['threads_jellyfish']
    conda:
        "../envs/jellyfish.yaml"
    shell:
        """
        set -e
        reads=$(cat {input.fof} | tr '\n' ' ')
        tmpfile=$(mktemp tmp/{wildcards.GE}_genomic_ccs_reads_XXXXXX.fa)
        zcat $reads > $tmpfile
        jellyfish count -C -m 21 -s 1000000000 -t {threads} $tmpfile -o {output.jf_k21}
        jellyfish histo -t {threads} {output.jf_k21} > {output.hist}
        """
        # # Read the FoF to get the paths to the ccs reads
        # reads=$(cat {input.fof} | tr '\n' ' ')
        
        # zcat $reads > tmp/{wildcards.GE}_genomic_ccs_reads.fa
        # jellyfish count -C -m 21 -s 1000000000 -t {threads} tmp/{wildcards.GE}_genomic_ccs_reads.fa -o {output.jf_k21}
        # jellyfish histo -t {threads} {output.jf_k21} > {output.hist}
        # rm tmp/{wildcards.GE}_genomic_ccs_reads.fa

# Rule genomescope
rule genomescope:
    input:
        rules.jellyfish_ccs.output,
        jf_k21 = "02_general/{GE}/03_genomescope/genomic_ccs_k21.jf",
        hist = "02_general/{GE}/03_genomescope/{GE}.hist"
    output:
        touch("02_general/{GE}/done__03_genomescope"),
        # Add other output files if needed
    params:
        script = "workflow/scripts/genomescope.R",
        outdir = "02_general/{GE}/03_genomescope",
        prefix = "{GE}_genomescope_"
    log:
        "logs/log__02_general_{GE}_03_genomescope"
    threads: config['threads_genomescope']
    conda:
        "../envs/r.yaml"
    shell:
        """
        Rscript {params.script} {input.hist} 21 10000 {params.outdir} 1000 1 > {log} 2>&1
        mv {params.outdir}/summary.txt {params.outdir}/{params.prefix}summary.txt
        mv {params.outdir}/model.txt {params.outdir}/{params.prefix}model.txt
        mv {params.outdir}/plot.png {params.outdir}/{params.prefix}plot.png
        mv {params.outdir}/plot.log.png {params.outdir}/{params.prefix}plot.log.png
        """
