rule repeatmasker:
    input:
        rules.rename_genome_assembly.output
    output:
        touch("04_annotation/01_repeatmasker/done__repeatmasker"),
        report("04_annotation/01_repeatmasker/" + ASSEMBLY_ID + ".fasta.tbl",
        caption = "../report/04_01_repeatmasker.rst", category = "04_annotation", subcategory = "01_repeatmasker")
    params:
        outdir = "04_annotation/01_repeatmasker",
        lib = config['repeatmasker_lib'],
        args = "-a -xsmall -no_is -s -e rmblast -dir ./ -gff"
    log:
        "logs/log__04_annotation_01_repeatmasker"
    threads: 
        config['threads_repeatmasker']
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        mkdir -p {params.outdir} && cd {params.outdir} && 
        RepeatMasker -pa {threads} {params.args} -lib ../../{params.lib} ../../{input} > ../../{log} 2>&1
        """


rule kraken2_silva_index:
    input:
    output:
        touch("04_annotation/done__04_kraken2_silva_index")
    threads: 1
    params:
        "--use-ftp --db tmp/silva --special silva"
    log:
        "logs/log__04_annotation_04_kraken2_silva_index"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2-build --use-ftp --db tmp/silva --download-taxonomy > {log} 2>&1 &&
        kraken2-build {params}  >> {log} 2>&1
        """

rule kraken2_silva:
    input:
        fasta = rules.rename_genome_assembly.output,
        silva = rules.kraken2_silva_index.output
    output:
        touch("04_annotation/done__kraken2_silva")
    params:
        args = "--db tmp/silva --confidence 0.5 --minimum-base-quality 20 --use-names --gzip-compressed",
        out = "04_annotation/04_kraken2_silva"
    threads:
        config['threads_kraken2']
    log:
        "logs/log__04_annotation_04_kraken2_silva"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 \
        {params.args} \
        --threads {threads} \
        --report {params.out}.report \
        {input.fasta} \
        --output {params.out} > {log} 2>&1
        """


###############################################################################

# MAKE CUSTOM KRAKEN DATABASE TO SEARCH FOR MITOCHONDRIAL READS

###############################################################################


rule kraken2_customDB:
    input:
    output:
        touch("04_annotation/done__06_kraken2_customDB")
    params:
        args = "--use-ftp --db kraken2_mitoDB",
        db = "workflow/resources/plan_mitogenomes_kraken_db.fa.gz"
    threads: 
        config['threads_kraken2']
    log:
        os.path.join(workflow.basedir, "logs/log__02_info_02_kraken2_index")
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        cd tmp/
        kraken2-build {params.args} --download-taxonomy > {log} 2>&1 && \
        kraken2-build {params.args} --download-library plasmid >> {log} 2>&1 && \
        zcat ../{params.db} > tmp.fa && \
        kraken2-build {params.args} --add-to-library tmp.fa >> {log} 2>&1 && \
        kraken2-build {params.args} --build --threads {threads} >> {log} 2>&1
        """
        

rule kraken2_mitosearch:
    input:
        k2 = rules.kraken2_customDB.output
    output:
        touch("04_annotation/done__06_kraken2_mitosearch")
    params:
        reads = " ".join(config["genome_ccs_files"]),
        args = "--db tmp/kraken2_mitoDB --confidence 0.5 --minimum-base-quality 20 --use-names --gzip-compressed",
        out = "04_annotation/06_kraken2/kraken2"
    threads:
        config['threads_kraken2']
    log:
        "logs/log__04_annotation_06_kraken2_mitosearch"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        cat {params.reads} > tmp/genomic_ccs_reads.fa.gz
        kraken2 {params.args} \
        --threads {threads} \
        --report {params.out}.report \
        tmp/genomic_ccs_reads.fa.gz \
        --output {params.out}  > {log} 2>&1
        """


###############################################################################

# BLAST SMALL LIST OF COMPLETE FLATWORM MITOCHONDRIAL GENOMES AGAINST ASSEMBLY

###############################################################################

rule make_genome_blastdb:
    input:
        rules.rename_genome_assembly.output
    output:
        touch("04_annotation/done__06_make_genome_blastdb")
    params:
        db_name = config['ucsc_id'],
        assembly_name = config['ucsc_id'] + ".fasta",
        db_prefix = config['ucsc_id'],
    log:
        os.path.join(workflow.basedir, "log/log__04_annotation_make_genome_blastdb")
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


rule blastn_mito_to_genome:
    input:
        rules.make_genome_blastdb.output,
        assembly = rules.rename_genome_assembly.output
    output:
        touch("04_annotation/done__06_blastn_mito"),
        tbl = report("04_annotation/06_blastn_mito/mito_to_" + config['ucsc_id'] + "blast.out",
        caption = "../report/04_06_blast_mito_to_genome.rst", category = "04_annotation", subcategory = "06_blastn_mito")
    params:
        db_name = "tmp/blastDB_genome/" + config['ucsc_id'],
        mito_db = "workflow/resources/plan_mitogenomes_kraken_db.fa.gz"
    threads:
        config['threads_blastn']
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        zcat {params.mito_db} > tmp/mito_query.fa
        blastn -db {params.db_name} \
        -num_threads {threads} \
        -query tmp/mito_query.fa \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qcovs sstart send evalue bitscore' > {output.tbl} &&
        rm tmp/mito_query.fa
        """

###############################################################################

# Find the coordinates of Ns

###############################################################################


rule find_Ns:
    input:
        rules.rename_genome_assembly.output
    output:
        touch("04_annotation/done__05_find_Ns"),
        bed = "04_annotation/05_ns/" + config['ucsc_id'] + "_N.bed"
    params:
        args = 'locate --bed -r -p "N+"',
    threads: 1
    log:
        "logs/log__04_annotation_05_find_Ns"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit {params.args} {input} > {output.bed}
        """
        

###############################################################################

# blast for telomeres | search for telomeres in repeatmasker output

###############################################################################


rule telomere_blast:
    input:
        rules.make_genome_blastdb.output
    output:
        touch("04_annotation/done__03_telomere_blast"),
        blastout = "04_annotation/03_telomere/" + config['ucsc_id'] + "_telomere.blastout",
        bed = "04_annotation/03_telomere/" + config['ucsc_id'] + "_telomere_blast.bed."
    params:
        db_name = "tmp/blastDB_genome/" + config['ucsc_id'],
        tel_db =  "workflow/resources/telomere_db.fa"
    threads: 1
    log:
        "logs/log__04_annotation_03_telomere_blast"
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        blastn -db {params.db_name} \
        -num_threads {threads} \
        -query {params.tel_db} -outfmt 6 > {output.blastout} &&
        cat {output.blastout} | sort -n -k2,2 -k9,9 | awk '{{printf "%s\\t%s\\t%s\\t%s\\t\\n", $2, $9, $10, $1}}' > {output.bed}
        """


rule telomere_in_repeatmasker:
    input:
        rules.repeatmasker.output
    output:
        touch("04_annotation/done__03_telomere_in_repeatmasker"),
        out = "04_annotation/03_telomere/" + config['ucsc_id'] + "_rm.out",
        gff = "04_annotation/03_telomere/" + config['ucsc_id'] + "_rm.gff"
    params:
        rm_out = "04_annotation/01_repeatmasker/" + config['ucsc_id'] + ".fasta.out",
        rm_gff = "04_annotation/01_repeatmasker/" + config['ucsc_id'] + ".fasta.out.gff"
    threads: 1
    log:
        "logs/log__04_annotation_03_telomere_in_repeatmasker"
    shell:
        """
        grep -e "TTAGGG" -e "SW" {params.rm_out} > {output.out}
        grep -e "TTAGGG" -e "gff" {params.rm_gff} > {output.gff}
        """


rule gem_index:
    input:
        rules.rename_genome_assembly.output
    output:
        touch("04_annotation/07_mappability/done__gem_index")
    params:
        outdir = "04_annotation/07_mappability",
        idx = "04_annotation/07_mappability/gem_index"
    log:
        "logs/log__04_annotation_07_mappability_gem_index"
    threads:
        config['threads_gem_index']
    shell:
        """
        export PATH=workflow/scripts/bin:$PATH
        mkdir -p {params.outdir} && 
        workflow/scripts/bin/gem-indexer -T {threads} -c dna -i {input} -o {params.idx} > {log} 2>&1
        """

rule gem:
    input:
        rules.rename_genome_assembly.output,
        rules.gem_index.output
    output:
        touch("04_annotation/07_mappability/done__gem")
    params:
        idx = "04_annotation/07_mappability/gem_index",
        out = "04_annotation/07_mappability/gem_out"
    log:
        "logs/log__04_annotation_07_mappability"
    threads:
        config['threads_gem']
    shell:
        """
        export PATH=workflow/scripts/bin:$PATH
        workflow/scripts/bin/gem-mappability -T {threads} -I {params.idx}.gem -l 100 -o {params.out} > {log} 2>&1
        """
 
rule gem_to_bw:
    input:
        gem = rules.gem.output,
        genome = rules.rename_genome_assembly.output
    output:
        touch("04_annotation/07_mappability/done__bw")
    params:
        idx = "04_annotation/07_mappability/gem_index",
        gem_out = rules.gem.params.out
    log:
        "logs/log__04_annotation_07_mappability_bw"
    shell:
        """
        export PATH=workflow/scripts/bin:$PATH
        samtools faidx {input.genome} && 
        cat {input.genome}.fai | awk '{{print$1"\t"$2}}' > {input.genome}.chrsizes &&
        workflow/scripts/bin/gem-2-wig -I {params.idx}.gem -i {params.gem_out}.mappability -o {params.gem_out} 2>> {log} && 
        cat {params.gem_out}.wig | awk -v OFS="\t" '{{print$1,$2,$5}}' > {params.gem_out}1.wig &&
        workflow/scripts/bin/wigToBigWig {params.gem_out}1.wig {input.genome}.chrsizes {params.gem_out}.bw > {log} 2>&1
        """



rule r04_annotation_done:
    input:
        rules.repeatmasker.output,
        rules.telomere_blast.output,
        rules.telomere_in_repeatmasker.output,
        rules.kraken2_silva.output,
        rules.find_Ns.output,
        rules.kraken2_mitosearch.output,
        rules.blastn_mito_to_genome.output,
        rules.gem_to_bw.output
    output:
        touch("04_annotation/done__04_annotation")
