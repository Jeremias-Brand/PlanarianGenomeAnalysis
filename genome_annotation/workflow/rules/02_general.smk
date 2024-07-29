C02 = "02 General"

rule merge_ccs:
    input:
        config['genome_ccs_files']
    output:
        config['merged_ccs_reads']
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule quast:
    input:
        rules.rename_genome_assembly.output
    output:
        "02_general/01_quast/report.tsv"
    log:
        "logs/log__02_general_01_quast"
    params:
        "02_general/01_quast"
    threads:
        config['threads_quast']
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast -o {params} -t {threads} {input} > {log} 2>&1
        """

rule kraken2_index:
    input:
    output:
        touch("02_general/02_kraken2/done__kraken2_index")
    params:
        "--use-ftp --db tmp/kraken2_index --download-library"
    threads: 1
    log:
        os.path.join(workflow.basedir, "logs/log__02_general_02_kraken2_index")
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        mkdir -p tmp/kraken2_index && cd tmp/kraken2_index &&
        wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20201202.tar.gz > {log} 2>&1 &&
        tar -xvf k2_standard_20201202.tar.gz >> {log} 2>&1
        """

rule kraken2:
    input:
        ccs = rules.merge_ccs.output,
        k2 = rules.kraken2_index.output
    output:
        touch("02_general/02_kraken2/done__kraken2")
    params:
        args = "--db tmp/kraken2_index --confidence 0.5 --minimum-base-quality 20 --use-names --gzip-compressed",
        out = "02_general/02_kraken2/kraken2"
    threads:
        config['threads_kraken2']
    log:
        "logs/log__02_general_02_kraken2"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 {params.args} \
        --threads {threads} \
        --report {params.out}.report \
        {input.ccs} \
        --output {params.out}  > {log} 2>&1
        """


###############################################################################

# GET NT BLAST DATABASE AND TAXONOMIC INFORMATION

###############################################################################


## TODO this is not checked with use-conda
rule get_nt:
    input:
    output:
        touch("02_general/done__02_nt_download")
    threads: config['threads_updateblast']
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        cd tmp/
        update_blastdb.pl --decompress --num_threads {threads} nt
        """

rule get_taxdump:
    output:
        "tmp/names.dmp",
        "tmp/nodes.dmp"
    threads: 1
    shell:
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P tmp/
        tar zxf tmp/taxdump.tar.gz -C tmp/ nodes.dmp names.dmp
        """

###############################################################################

# CREATE DIAMOND DATABASE FROM UNIPROT

###############################################################################


rule get_uniprot:
    output:
        touch("02_general/done__get_uniprot")
    params: 
        version = config["uniprot_version"]
    threads: 1
    shell:
        """
        mkdir -p tmp/uniprot
        cd tmp/uniprot && 
        wget -c {params.version} && 
        tar zxvf *.tar.gz
        """

rule prepare_uniprot:
    input:
        "02_general/done__get_uniprot"
    output:
        touch("02_general/done__04_prepare_uniprot"),
        "tmp/uniprot/uniprot_ref_proteomes.fasta",
        "tmp/uniprot/uniprot_ref_proteomes.taxids"
    log:
        "logs/prepare_uniprot.log"
    threads: config['threads_prepare_uniprot']
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        cd tmp/uniprot/
        # Unpack protein FASTAs for each kingdom
        find . -name *.fasta.gz | grep -v 'DNA' | grep -v 'additional' | sed -E 's/(.+)/gunzip \\1/' > gunzip.cmd
        parallel -j{threads} -a gunzip.cmd
        find . -name *.idmapping.gz |  sed -E 's/(.+)/gunzip \\1/' > gunzip_id.cmd
        parallel -j{threads} -a gunzip_id.cmd
        # Concatenate all protein sequences into uniprot_ref_proteomes.fasta
        cat */*/*.fasta > uniprot_ref_proteomes.fasta
        # Subset mapping file to only contain NCBI TaxID entries
        cat */*/*.idmapping | grep "NCBI_TaxID" > uniprot_ref_proteomes.taxids
        # Simplify sequence IDs
        cat uniprot_ref_proteomes.fasta | sed -r 's/(^>sp\|)|(^>tr\|)/>/g' | cut -f1 -d"|" > temp; mv temp uniprot_ref_proteomes.fasta
        """

rule makedb_diamond:
    input:
        "02_general/done__04_prepare_uniprot"
    output:
        touch("02_general/done__makedb_diamond"),
    log:
        "logs/02__makedb_diamond.log"
    threads: 128
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        cd tmp/uniprot/
        diamond makedb --in uniprot_ref_proteomes.fasta -d uniprot_ref_proteomes
        """


###############################################################################

# SEARCH NT AND UNIPROT WITH THE ASSEMBLY

###############################################################################


rule blastn_assembly:
    input:
        "02_general/done__02_nt_download",
        assembly = rules.rename_genome_assembly.output
    output:
        touch("02_general/done__blastn"),
        blastout = "02_general/" + config["ucsc_id"] + "_ccs.blastout"
    threads: config['threads_blastn']
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        blastn \
        -query {input.assembly} \
        -db tmp/nt \
        -outfmt '6 qseqid staxids bitscore std' \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads {threads} > {output.blastout}
        """


rule diamond_assembly:
    input:
        rules.makedb_diamond.output,
        assembly = rules.rename_genome_assembly.output
    output:
        "02_general/" + config['ucsc_id']+ "_diamond.out"
    threads: config['threads_diamond']
    params:
        "tmp/uniprot/uniprot_ref_proteomes.dmnd"
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        diamond blastx \
        --query {input.assembly} \
        --db {params} \
        --outfmt 6 \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads {threads} \
        --out {output}
        """


###############################################################################

# CREATE BAM FILE WITH GENOMIC READS

###############################################################################


rule map_ccs_genome:
    input:
        rules.rename_genome_assembly.output
    output:
        touch("02_general/done__map_ccs_genome"),
        bam = "02_general/" + config["ucsc_id"] + "_ccs.bam"
    params:
        ID = config["ucsc_id"],
        sample = "genomic_CCS",
        reads = " ".join(config["genome_ccs_files"])
    log:
        "logs/log__02_map_ccs_genome"
    threads: config['threads_pbmm2']
    conda:
        "../envs/pbmm2.yaml"
    shell:
        """
        zcat {params.reads} > tmp/genomic_ccs_reads.fa
        pbmm2 align {input} tmp/genomic_ccs_reads.fa {output.bam} --preset CCS --sort --rg  '@RG\tID:{params.ID}\tSM:{params.sample}' > {log} 2>&1
        rm tmp/genomic_ccs_reads.fa
        """


###############################################################################

# GET BLOBTOOLS AND PREPARE EVERYTHING FOR PLOTTING

###############################################################################


rule get_blobtools:
    input:
    output:
        "tmp/blobtools/blobtools"
    params:
        blob_link = "https://github.com/DRL/blobtools/archive/refs/tags/blobtools_v1.1.1.tar.gz",
        blob_name = "blobtools_v1.1.1.tar.gz"
    threads: 1
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        wget {params.blob_link} -P tmp/ 
        mkdir -p tmp/b_tmp
        tar zxf tmp/{params.blob_name} -C tmp/b_tmp
        mkdir -p tmp/blobtools
        mv tmp/b_tmp/*/* tmp/blobtools/
        """


rule taxify_diamond:
    input:
        blobtools = "tmp/blobtools/blobtools",
        diamond_out = "02_general/" + config['ucsc_id'] + "_diamond.out"
    output:
        diamond_tax = "02_general/" + config['ucsc_id'] + "_diamond.taxified.out"
    threads: 1
    log:
        "logs/02__taxify_diamond.log"
    params:
        pre = "02_general/",
        tax = "tmp/uniprot/uniprot_ref_proteomes.taxids"
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        {input.blobtools} taxify -f {input.diamond_out} -m {params.tax} -s 0 -t 1 -o {params.pre} 2> {log}
        """


rule blobtools_map2cov:
    input:
        blobtools = "tmp/blobtools/blobtools",
        assembly = rules.rename_genome_assembly.output,
        bam = "02_general/" + config['ucsc_id'] + "_ccs.bam"
    output:
        "02_general/" + config['ucsc_id'] + "_ccs.bam.cov"
    params:
        prefix = "02_general/"
    log:
        "logs/02__blobtools_map2cov.log"
    conda:
        "../envs/blobtools.yaml"
    threads: 1
    shell:
        """
        {input.blobtools} map2cov -i {input.assembly} -b {input.bam} -o {params.prefix} 2> {log}
        """


rule blobtools_create:
    input:
        blobtools = "tmp/blobtools/blobtools",
        assembly = rules.rename_genome_assembly.output,
        cov = "02_general/" + config['ucsc_id'] + "_ccs.bam.cov",
        blast_out = "02_general/" + config['ucsc_id'] + "_ccs.blastout",
        diamond_tax = "02_general/" + config['ucsc_id'] + "_diamond.taxified.out",
        names = "tmp/names.dmp",
        nodes = "tmp/nodes.dmp"
    output:
        db = "02_general/blobplot/" + config['ucsc_id'] + ".blobDB.json"
    log:
        "logs/02__blobtools_create.log"
    params:
        prefix = "02_general/blobplot/" + config['ucsc_id']
    conda:
        "../envs/blobtools.yaml"
    threads: 1
    shell:
        """
        {input.blobtools} create \
        -i {input.assembly} \
        -c {input.cov} \
        -t {input.blast_out} \
        -t {input.diamond_tax} \
        -o {params.prefix} \
        --names {input.names} \
        --nodes {input.nodes} &> {log}
        """

###############################################################################

# MAKE BLOBPLOT

###############################################################################


rule blobplot:
    input:
        blobtools = "tmp/blobtools/blobtools",
        db = rules.blobtools_create.output
    output:
        stats = "02_general/blobplot/" + config['ucsc_id'] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.stats.txt",
        p1  = "02_general/blobplot/" + config['ucsc_id'] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.cov0.png",
        p2 = "02_general/blobplot/" + config['ucsc_id'] + ".blobDB.json.bestsum.phylum.p8.span.100.blobplot.read_cov.cov0.png",
        done = touch("02_general/done__02_blobplot")
    params:
        prefix = "02_general/blobplot/"
    log:
        "logs/02__02_blobplot.log"
    threads: 1
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        echo '##### VIEW #####' > {log}
        {input.blobtools} view \
        -i {input.db} \
        -o {params.prefix} &>> {log}
        echo '##### PLOT #####' >> {log}
        
        {input.blobtools} plot \
        -i {input.db} \
        -o {params.prefix} &>> {log}
        """

###############################################################################

# GET K-MER SPECTRUM AND RUN GENOMESCOPE

###############################################################################


rule jellyfish_ccs:
    input:
        config['genome_ccs_files']
    output:
        touch("02_general/done__03_jellyfish"),
        jf_k21 = "02_general/03_genomescope/genomic_ccs_k21.jf",
        hist = "02_general/03_genomescope/" + config['ucsc_id']+ ".hist"
    params:
        reads = " ".join(config["genome_ccs_files"])
    threads: 
        config['threads_jellyfish']
    conda:
        "../envs/jellyfish.yaml"
    shell:
        """
        zcat {params.reads} > tmp/genomic_ccs_reads.fa
        jellyfish count -C -m 21 -s 1000000000 -t {threads} tmp/genomic_ccs_reads.fa -o {output.jf_k21}
        jellyfish histo -t {threads} {output.jf_k21} > {output.hist}
        rm tmp/genomic_ccs_reads.fa
        """


rule genomescope:
    input:
        rules.jellyfish_ccs.output,
        jf_k21 = "02_general/03_genomescope/genomic_ccs_k21.jf",
        hist = "02_general/03_genomescope/" + config['ucsc_id']+ ".hist"
    output:
        touch("02_general/done__03_genomescope"),
        
        #report([ "02_general/03_genomescope/" + config['ucsc_id'] + "_genomescope_" +  i for i in ["summary.txt", "model.txt", "plot.png", "plot.log.png"]],
        report([ "02_general/03_genomescope/%s_genomescope_%s" % (config['ucsc_id'], i) 
        for i in ["summary.txt", "model.txt", "plot.png", "plot.log.png"]],
        caption = "../report/02_genomescope.rst", category="02 General", subcategory="genomescope")
    params:
        script = "workflow/scripts/genomescope.R",
        outdir = "02_general/03_genomescope",
        prefix = config['ucsc_id'] + "_genomescope_"
    log:
        "logs/log__02_general_03_genomescope"
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

rule r02_general_done:
    input:
 #       rules.merge_ccs.output,
 #       rules.quast.output,
 #       rules.kraken2.output,
        #rules.blastn_assembly.output,
        #rules.diamond_assembly.output,
#        rules.map_ccs_genome.output,
#        rules.blobplot.output,
        rules.genomescope.output
    output:
        touch("02_general/done__02_general")
 
