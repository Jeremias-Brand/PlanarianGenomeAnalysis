# some of this is blanked out because we are running this pipeline as a subworkflow
#configfile: "config/config.yaml"

onsuccess:
    print("Workflow finished, no error!"),
    shell("mail -s 'snakemake finished OK' jbrand@mpibpc.mpg.de < {log}")

onerror:
    print("An error occurred"),
    shell("mail -s 'snakemake finished NOT OK' jbrand@mpibpc.mpg.de < {log}")

ref_genome = config["ref_genome"]
scaffolds = []
# collect the contig names
with open(ref_genome, "rt") as infile:
    for line in infile:
        line = line.strip()
        if line.startswith(">"):
            line = line.split(" ")[0]
            scaffolds.append(line[1:])

# To make the location of the pipeline more flexible we give a base directory here

vardir = ""

SEdir = vardir + "SE/"
PEdir = vardir + "PE/"
GATKdir = vardir + "PE/GATK_" + ref_genome[4:-3] + "/"
ref_name = ref_genome[4:-3]
LOGdir = os.path.join(workflow.basedir, "logs/log__PE_GATK_" + ref_genome + "_")


rule all_PE_GATK:
    input:
        GATKdir + "done__PE_GATK_joint"
    output:
        touch("06_variation/done__PE_GATK_all")


rule reset_genotyping:
    params:
        SE = SEdir,
        PE = PEdir
    shell:
        """
        rm -f 06_variation/done__01_all
        rm -rdf {params.SE}QC
        rm -rdf {params.SE}{{02,03,04,05,06}}*
        rm -rdf {params.PE}QC
        rm -rdf {params.PE}{{02,03,04,05,06}}*
        """

##############################################################
#   INDEX THE REFERENCE GENOME
############################################################## 

rule PE_GATK_index_ref_genome:
    input:
        ref_genome
    output:
    #TODO set this to temp()
        ref_genome + ".bwt"
    conda:
        "../envs/06_bwa.yaml"
    log:
        LOGdir + "index_ref_genome_" + ref_genome
    shell:
        """
        bwa index -a bwtsw {input} &> {log}
        """


rule PE_GATK_samdict_ref_genome:
    input:
        ref_genome
    output:
    #TODO set this to temp()
        ref_genome + ".dict"
    conda:
        "../envs/06_samtools.yaml"
    log:
        LOGdir + "samdict_ref_genome_" + ref_genome
    shell:
        """
        samtools dict {input} > {output} 2> {log}
        """


rule PE_GATK_seqdict_ref_genome:
    input:
        ref_genome
    output:
    #TODO set this to temp()
        ref_genome + ".seqdic.done"
    conda:
        "../envs/gatk.yaml"
    log:
        LOGdir + "seqdict_ref_genome_" + ref_genome
    shell:
        """
        gatk CreateSequenceDictionary -R {input} 2> {log} && touch {output}
        """


rule PE_GATK_faidx_ref_genome:
    input:
        ref_genome
    output:
        ref_genome + ".fai"
    conda:
        "../envs/06_samtools.yaml"
    log:
        LOGdir + "faidx_ref_genome_" + ref_genome
    shell:
        """
        samtools faidx {input} > {output} 2> {log}
        """

###############################################################################

# fastq2ubam PE

###############################################################################

rule PE_GATK_fastq2ubam:
    input:
        r1 = PEdir + "02_trim/{PE}_1.fastq.gz",
        r2 = PEdir + "02_trim/{PE}_2.fastq.gz",
    output:
        temp(GATKdir + "03_bam/{GE}__{PE}_unaligned_reads.bam")
    conda:
        "../envs/06_picard.yaml"
    threads: 8
    log:
        LOGdir + "PE_GATK_fastq2ubam_{GE}__{PE}"
    shell:
        """
        picard FastqToSam \
         F1={input.r1} \
         F2={input.r2} \
         OUTPUT={output} \
         PLATFORM=illumina \
         SAMPLE_NAME={wildcards.PE} &> {log}
        """


###############################################################################

# ubam2gvcf PE

###############################################################################

rule PE_GATK_SamToFastq_bwa:
# interleave for paired!
    input:
        ubam = rules.PE_GATK_fastq2ubam.output,
        ref_genome = ref_genome,
        ref_genome_index = ref_genome + ".bwt"
    output:
        temp(GATKdir + "03_bam/{GE}__{PE}_unmerged.bam")
    conda:
        "../envs/06_picard.yaml"
    log:
        l1 = LOGdir + "PE_GATK_SamToFastq_bwa_l1_{GE}__{PE}",
        l2 = LOGdir + "PE_GATK_SamToFastq_bwa_l2_{GE}__{PE}",
        l3 = LOGdir + "PE_GATK_SamToFastq_bwa_l3_{GE}__{PE}"
    threads: config["bwa_threads"]
    shell:
        """
        picard SamToFastq \
        INPUT={input.ubam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true 2> {log.l1}| \
        bwa mem -K 100000000 -p -v 3 -t {threads} -Y {input.ref_genome} /dev/stdin 2> {log.l2} | \
        samtools view -1 - > {output} && [[ -s {output} ]]
        """

rule PE_GATK_MergeBamAlignment:
    # This adds alignment information stored in the unaligned with the aligned bam file.
    # https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
    input:
        unmerged_bam = rules.PE_GATK_SamToFastq_bwa.output,
        ubam = GATKdir + "03_bam/{GE}__{PE}_unaligned_reads.bam",
        ref_genome = ref_genome, 
        ref_genome_sam_index = ref_genome + ".dict"
    output:
        temp(GATKdir + "03_bam/{GE}__{PE}_aligned_unsorted.bam")
    conda:
        "../envs/06_gatk.yaml"
    log:
        LOGdir + "PE_GATK_MergeBamAlignment_{GE}__{PE}"
    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" \
        MergeBamAlignment \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ALIGNED_BAM {input.unmerged_bam}  \
        --CREATE_INDEX true \
        --UNMAPPED_BAM {input.ubam} \
        --OUTPUT {output} \
        --REFERENCE_SEQUENCE {input.ref_genome} \
        --PAIRED_RUN false \
        --SORT_ORDER "unsorted" \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --CLIP_ADAPTERS false \
        --CLIP_OVERLAPPING_READS true \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --UNMAP_CONTAMINANT_READS true &> {log}
        """


rule PE_GATK_MarkDuplicates:
    input:
        rules.PE_GATK_MergeBamAlignment.output
    output:
        bam = temp(GATKdir + "03_bam/{GE}__{PE}_aligned_unsorted_duplicates_marked.bam"),
        metrics_filename = GATKdir + "QC/{GE}__{PE}.dup_metrics"
    conda:
        "../envs/06_gatk.yaml"
    log:
        LOGdir + "PE_GATK_MarkDuplicates_{GE}__{PE}"
    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m -XX:+UseParallelGC -XX:ParallelGCThreads=2" \
        MarkDuplicates \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics_filename} \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        --ASSUME_SORT_ORDER "queryname" \
        --CREATE_MD5_FILE true &> {log}
        """

rule PE_GATK_SortAndFixTags:
    input:
        bam = rules.PE_GATK_MarkDuplicates.output.bam,
        ref_genome = ref_genome
    output:
        GATKdir + "03_bam/{GE}__{PE}_aligned_duplicates_marked_sorted.bam"
    conda:
        "../envs/06_gatk.yaml"
    log:
        LOGdir + "PE_GATK_SortAndFixTags_{GE}__{PE}"
    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" \
        SortSam \
        --INPUT {input.bam} \
        --OUTPUT /dev/stdout \
        --SORT_ORDER "coordinate" \
        --CREATE_INDEX false \
        --CREATE_MD5_FILE false 2> {log}\
        | \
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms500m" \
        SetNmMdAndUqTags \
        --INPUT /dev/stdin \
        --OUTPUT {output} \
        --CREATE_INDEX true \
        --CREATE_MD5_FILE true \
        --REFERENCE_SEQUENCE {input.ref_genome} &>> {log}
        """

rule PE_GATK_HaplotypeCaller_gvcf:
    input:
        bam = rules.PE_GATK_SortAndFixTags.output,
        ref_genome = ref_genome,
        ref_genome_index = ref_genome + ".fai",
        gatk_dict = ref_genome + ".seqdic.done",
        sam_dict = ref_genome + ".dict"
    output:
        GATKdir + "04_gvcf/{GE}__{PE}.g.vcf.gz"
    conda:
        "../envs/06_gatk.yaml"
    log:
        LOGdir + "PE_GATK_HaplotypeCaller_gvcf_{GE}__{PE}"
    threads: config["threads_HaplotypeCaller"]
    shell:
        """
        gatk --java-options "-Xmx20G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP \
        --native-pair-hmm-threads {threads} \
        -R {input.ref_genome} \
        -I {input.bam} \
        -O {output} \
        -contamination 0 -ERC GVCF &> {log}
        """

rule PE_GATK_HaplotypeCaller_vcf:
    input:
        bam = rules.PE_GATK_SortAndFixTags.output,
        ref_genome = ref_genome,
        ref_genome_index = ref_genome + ".fai",
        gatk_dict = ref_genome + ".seqdic.done",
        sam_dict =  ref_genome + ".dict"
    output:
        GATKdir + "06_vcf/single/{GE}__{PE}.vcf.gz"
    conda:
        "../envs/06_gatk.yaml"
    log:
        LOGdir + "PE_GATK_HaplotypeCaller_vcf_{GE}__{PE}"
    threads: config["threads_HaplotypeCaller"]
    shell:
        """
        gatk --java-options "-Xmx20G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP \
        --native-pair-hmm-threads {threads} \
        -R {input.ref_genome} \
        -I {input.bam} \
        -O {output} \
        -contamination 0 &> {log}⁄
        """

rule PE_GATK_Gather_single_vcfs:
    input:
        expand(GATKdir + "06_vcf/single/" + ref_genome + "__{PE}.vcf.gz", PE = config["PE"])
    output:
        file = GATKdir + "06_vcf/single_all_raw.vcf.gz",
        index = GATKdir + "06_vcf/single_all_raw.vcf.gz.tbi"
    conda:
        "../envs/gatk.yaml"
    params:
        " -I ".join(GATKdir + "06_vcf/single/" + s + ".vcf.gz" for s in config["PE"])
    log:
        out = LOGdir + "PE_GATK_Gather_single_vcfs",
        index = LOGdir + "PE_GATK_Gather_single_vcfs_index"
    shell:
        """
        gatk --java-options "-Xmx24g -Xms24g" \
        MergeVcfs \
        -I {params} \
        -O {output.file} &> {log.out}

        gatk --java-options "-Xmx24g -Xms24g" \
        IndexFeatureFile \
        -I {output.file} \
        -O {output.index} &> {log.index}
        """

# We have made SNP calls with each individual called .vcf files
# now we also potentially do joint genotyping. 
# for this we need to gather and jointly genotype all the gvcf files
GATKdir = vardir + "PE/GATK_" + ref_genome[4:-3] + "/"


###############################################################################

# genotyping PE

###############################################################################

rule PE_GATK_WriteSampleMap:
    input:
        #"metadata.tsv",
        expand(GATKdir + "04_gvcf/" + ref_name + "__{PE}.g.vcf.gz", PE=config["PE"])
    output:
        GATKdir + "PE_GATK_gvcf_sample_map.txt"
    run:
        #import pandas as pd
        #samples = pd.read_csv("metadata.tsv",sep='\t').set_index(["prefix"], drop=False)
        files = [GATKdir + "04_gvcf/" + ref_name + "__" + i + ".g.vcf.gz" for i in config["PE"]]
        sample = [os.path.basename(i).split(".")[0] for i in files]
        out = open(output[0],'w')
        for c1,c2 in zip(sample, files):
            print("%s\t%s"%(c1,c2),file = out)


rule PE_GATK_GenomicsDBImport:
    input:
        GATKdir + "PE_GATK_gvcf_sample_map.txt"
    output:
        db = directory(GATKdir + "05_genomicsDB/{scaffold}.db"),
        tar = GATKdir + "05_genomicsDB/{scaffold}.db.tar"
    conda:
        "../envs/gatk.yaml"
    log:
        db = LOGdir + "PE_GATK_GenomicsDBImport_{scaffold}",
        tar = LOGdir + "PE_GATK_GenomicsDBImport_tar_{scaffold}"
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output.db} \
        --L {wildcards.scaffold} \
        --sample-name-map {input} \
        --reader-threads 5 \
        --batch-size 50 &> {log.db}

        tar -cf {output.tar} {output.db} &> {log.tar}
        """

rule PE_GATK_GenotypeGVCFs:
    input:
        db = GATKdir + "05_genomicsDB/{scaffold}.db",
        ref_genome = ref_genome
    output:
        GATKdir + "06_vcf/scaffolds/{scaffold}.g.vcf.gz"
    conda:
        "../envs/gatk.yaml"
    log:
        LOGdir + "PE_GATK_GenotypeGVCFs_{scaffold}"
    shell:
        """
        gatk --java-options "-Xmx8g -Xms8g" \
        GenotypeGVCFs \
        -R {input.ref_genome} \
        -O {output} \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        -V gendb://{input.db} \
        -L {wildcards.scaffold} &> {log}
        """

rule PE_GATK_Gather_joint_vcfs:
    input:
        expand(GATKdir + "06_vcf/scaffolds/{scaffold}.g.vcf.gz",scaffold = scaffolds)
    output:
        file = GATKdir + "06_vcf/joint_all_raw.vcf.gz",
        index = GATKdir + "06_vcf/joint_all_raw.vcf.gz.tbi"
    conda:
        "../envs/gatk.yaml"
    params:
        " -I ".join(GATKdir + "06_vcf/scaffolds/" + s + ".g.vcf.gz" for s in scaffolds)
    log:
        out = LOGdir + "PE_GATK_GatherVcfs",
        index = LOGdir + "PE_GATK_GatherVcfs_index"
    shell:
        """
        gatk --java-options "-Xmx24g -Xms24g" \
        MergeVcfs \
        -I {params} \
        -O {output.file} &> {log.out}

        gatk --java-options "-Xmx24g -Xms24g" \
        IndexFeatureFile \
        -I {output.file} \
        -O {output.index} &> {log.index}
        """

###############################################################################

# filtering PE

###############################################################################

rule PE_GATK_single_selectSNPs:
    input:
        GATKdir + "06_vcf/single_all_raw.vcf.gz"
    output:
        GATKdir + "06_vcf/single_all_raw_snp.vcf.gz"
    conda:
        "../envs/06_gatk.yaml"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type SNP \
        -O {output}
        """


rule PE_GATK_joint_selectSNPs:
    input:
        GATKdir + "06_vcf/joint_all_raw.vcf.gz"
    output:
        GATKdir + "06_vcf/joint_all_raw_snp.vcf.gz"
    conda:
        "../envs/06_gatk.yaml"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type SNP \
        -O {output}
        """


rule PE_GATK_single_selectINDELs:
    input:
        GATKdir + "06_vcf/single_all_raw.vcf.gz"
    output:
        GATKdir + "06_vcf/single_all_raw_indel.vcf.gz"
    conda:
        "../envs/06_gatk.yaml"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type INDEL \
        -O {output}
        """


rule PE_GATK_joint_selectINDELs:
    input:
        GATKdir + "06_vcf/joint_all_raw.vcf.gz"
    output:
        GATKdir + "06_vcf/joint_all_raw_indel.vcf.gz"
    conda:
        "../envs/06_gatk.yaml"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type INDEL \
        -O {output}
        """


rule PE_GATK_single_gather:
    input:
        rules.PE_GATK_single_selectSNPs.output,
        rules.PE_GATK_single_selectINDELs.output
    output:
        touch(GATKdir + "done__PE_GATK_single")


rule PE_GATK_joint_gather:
    input:
        rules.PE_GATK_joint_selectSNPs.output,
        rules.PE_GATK_joint_selectINDELs.output
    output:
        touch(GATKdir + "done__PE_GATK_joint")


