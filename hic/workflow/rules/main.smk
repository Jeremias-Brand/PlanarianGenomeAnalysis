# generate target list 
target_dict = {}

with open(config["combinations"], "r") as infile:
    for line in infile:
        key, value = line.strip().split("\t")
        if key in target_dict:
            target_dict[key].append(value)
        else:
            target_dict[key] = [value]


#########################################################################################


rule map_mates:
    input:
        genome = config["genome_dir"] + '{GE}.fa',
        R1 = config["omnic_dir"] + '{HIC}_R1.fq.gz',
        R2 = config["omnic_dir"] + '{HIC}_R2.fq.gz'
    output:
        touch("01_bam/done__map_mates_{GE}__{HIC}"),
        sam = temp('01_bam/{GE}__{HIC}_aligned.sam')
    params:
        dir = os.path.join(workflow.basedir, config["omnic_dir"]),
        opt = config['opt_bwa']
    log:
        "01_bam/log__map_mates_{GE}__{HIC}"
    threads: config['threads_bwa']
    conda:
        "../envs/omnic.yaml"
    shell:
        """
        if [ ! -f {input.genome}.pac ]; then
            bwa index {input.genome} 2> {log}
        else
            echo "bwa index exists" > {log}
        fi
        bwa mem {params.opt} -t {threads} {input.genome} {input.R1} {input.R2} -o {output.sam} 2>> {log}
        """

rule make_chromlen:
    input:
        index = config["genome_dir"] + '{GE}.fa.fai'
    output:
        config["genome_dir"] + '{GE}.chromsize'
    shell:
        """
        cut -f1,2 {input} > {output}
        """

rule pairtools_parse:
    input:
        rules.map_mates.output,
        index = config["genome_dir"] + '{GE}.fa.fai',
        sam = rules.map_mates.output.sam
    output:
        touch("01_bam/done__pairtools_parse_{GE}__{HIC}"),
        out = temp("01_bam/{GE}__{HIC}.pairsam")
    params:
        opt = config['opt_pairtools_parse'],
    log:
        "01_bam/log__pairtools_parse_{GE}__{HIC}"
    threads: config['threads_bwa']
    conda:
        "../envs/omnic.yaml"
    shell:
        """
        pairtools parse {params.opt} --nproc-in {threads} --nproc-out {threads} --chroms-path {input.index} {input.sam} > {output.out} 2> {log}
        """


rule pairtools_sort:
    input:
        rules.pairtools_parse.output.out
    output:
        touch("01_bam/done__pairtools_sort_{GE}__{HIC}"),
        out = temp("01_bam/{GE}__{HIC}_sorted.pairsam")
    params:
        tmp = "tmp/"
    log:
        "01_bam/log__pairtools_sort_{GE}__{HIC}"
    threads: config['threads_bwa']
    conda:
        "../envs/omnic.yaml"
    shell:
        """
        pairtools sort --nproc {threads} --tmpdir={params.tmp} {input} > {output.out} 2> {log}
        """


rule pairtools_dedup:
    input:
        rules.pairtools_sort.output.out
    output:
        touch("01_bam/done__pairtools_dedup_{GE}__{HIC}"),
        out = temp("01_bam/{GE}__{HIC}_sorted_dedup.pairsam"),
        stats =  "01_bam/{GE}__{HIC}_stats.txt"
    log:
        "01_bam/log__pairtools_dedup_{GE}__{HIC}"
    threads: config['threads_bwa']
    conda:
        "../envs/omnic.yaml"
    shell:
        """
        pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups --output-stats {output.stats} \
        --output {output.out} {input} 2> {log}
        """


rule pairtools_split:
    input:
        rules.pairtools_dedup.output.out
    output:
        touch("01_bam/done__pairtools_split_{GE}__{HIC}"),
        pairs = "01_bam/{GE}__{HIC}_mapped.pairs",
        bam =  temp("01_bam/{GE}__{HIC}_unsorted.bam")
    log:
        "01_bam/log__pairtools_dedup_{GE}__{HIC}"
    threads: config['threads_bwa']
    conda:
        "../envs/omnic.yaml"
    shell:
        """
        pairtools split --nproc-in {threads} --nproc-out {threads} \
        --output-pairs {output.pairs} \
        --output-sam {output.bam} {input}
        """


rule samtools_sort_index:
    input:
        bam = rules.pairtools_split.output.bam
    output:
        touch("01_bam/done__samtools_sort_index_{GE}__{HIC}"),
        bam =  "01_bam/{GE}__{HIC}_mapped_sorted.bam"
    params:
        tmp = "tmp/"
    log:
        "01_bam/log__samtools_sort_index_{GE}__{HIC}"
    threads: config['threads_bwa']
    conda:
        "../envs/omnic.yaml"
    shell:
        """
        samtools sort -@{threads} -T {params.tmp}/tempfile_$( date +%s%N ).bam -o {output.bam} {input.bam} 2> {log} && \
        samtools index {output.bam} 2> {log}
        """


rule pairs_qc:
    input:
        stats =  "01_bam/{GE}__{HIC}_stats.txt",
        bam = rules.samtools_sort_index.output.bam
    output:
        touch("02_qc/done__pairs_qc_{GE}__{HIC}"),
        out = "02_qc/{GE}__{HIC}_stats.txt",
        preseq = "02_qc/{GE}__{HIC}_preseq.txt"
    params:
        script = os.path.join(workflow.basedir, 'workflow/scripts/get_qc.py')
    log:
        "02_qc/log__pairs_qc_{GE}__{HIC}"
    threads: config['threads_bwa']
    conda:
        "../envs/omnic.yaml"
    shell:
        """
        preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output {output.preseq} {input.bam}
        python3 {params.script} -p {input.stats} -d {output.preseq} > {output.out} 2> {log} 
        """

# infer target list based on the combination file (see at the beginning of the script)
qc_list = []
for key in target_dict.keys():
    #print("Mapping these data: " + str(target_dict[key]) + " to this reference: " + key)
    for val in target_dict[key]:
        qc_list.append("02_qc/done__pairs_qc_" + str(key) +"__"+ str(val))

rule all_pairs:
    input:
        qc_list
    output:
        touch("02_qc/done__all_pairs")


rule pair_index:
    input:
        rules.pairtools_split.output.pairs
    output:
        touch("03_pairs/done__pairs_index_{GE}__{HIC}"),
        pair = "03_pairs/{GE}__{HIC}_mapped.pairs.gz"
    log:
        "03_pairs/log__pair_index_{GE}__{HIC}"
    threads: config['threads_bgzip']
    conda:
        "../envs/cooler.yaml"
    shell:
        """
        bgzip -@ {threads} -c {input} > {output.pair} && \
        pairix {output.pair}
        """

rule cooler_maxres:
    input:
        pair = rules.pair_index.output.pair,
        chromsize = rules.make_chromlen.output
    output:
        touch("04_cool/done__cooler_maxres_{GE}__{HIC}"),
        cool = "04_cool/{GE}__{HIC}_maxres_" + config["maxres_cool"] + ".cool"
    params:
        maxres = config["maxres_cool"]
    log:
        "04_cool/log__cooler_maxres_{GE}__{HIC}"
    threads: config['threads_cool']
    conda:
        "../envs/cooler.yaml"
    shell:
        """
        cooler cload pairix -p {threads} {input.chromsize}:{params.maxres} {input.pair} {output.cool} 2> {log}
        """


rule cooler_exres:
    input:
        pair = rules.pair_index.output.pair,
        chromsize = rules.make_chromlen.output
    output:
        touch("04_cool/done__cooler_exres_{GE}__{HIC}"),
        cool = "04_cool/{GE}__{HIC}_exres_" + config["exres_cool"] + ".cool"
    params:
        exres = config["exres_cool"]
    log:
        "04_cool/log__cooler_exres_{GE}__{HIC}"
    threads: config['threads_cool']
    conda:
        "../envs/cooler.yaml"
    shell:
        """
        cooler cload pairix -p {threads} {input.chromsize}:{params.exres} {input.pair} {output.cool} 2> {log}
        """


rule cooler_balance:
    # This balances an existing .cool file
    input:
        cool = "04_cool/{GE}__{HIC}_exres_" + config["exres_cool"] + ".cool"
    output:
        touch("04_cool/done__cooler_balance_{GE}__{HIC}")
        #cool = "04_cool/{GE}__{HIC}_exres_" + config["exres_cool"] + "_balanced.cool"
    params:
        exres = config["exres_cool"]
    log:
        "04_cool/log__cooler_balance_{GE}__{HIC}"
    threads: config['threads_cool']
    conda:
        "../envs/cooler.yaml"
    shell:
        """
        cooler balance -p {threads} {input.cool} 2> {log}
        """


rule cooler_zoomify:
    input:
        rules.cooler_maxres.output.cool
    output:
        touch("04_cool/done__cooler_zoomify_{GE}__{HIC}"),
        mcool = "04_cool/{GE}__{HIC}.mcool"
    params:
        opt = config["zoomify_opt"]
    log:
        "04_cool/log__cooler_zoomify_{GE}__{HIC}"
    threads: config['threads_cool']
    conda:
        "../envs/cooler.yaml"
    shell:
        """
        cooler zoomify -p {threads} {params.opt} -o {output.mcool} {input} 2> {log}
        """


# infer target list based on the combination file (see at the beginning of the script)
# TODO work on the expectation matrix at the moment it is turned off.
cool_list = []
for key in target_dict.keys():
    #print("Generating HiC cool matrix for this data : " + str(target_dict[key]) + " on this reference:  " + key)
    for val in target_dict[key]:
        cool_list.append("04_cool/" + str(key) +"__"+ str(val) + ".mcool")
        #cool_list.append("04_cool/done__cooltools_exp_" + str(key) +"__"+ str(val))

rule all_cool:
    input:
        cool_list
    output:
        touch("04_cool/done__all_cool")


# GENERATE EXPECTATION FILES FROM THE COOL FILES WITH COOLTOOLS
# bedtools.makewindows.-g.lib/schMedA2.chromsize -w 1000000 > lib/schMedA2.1Mb.bed
# bedtools makewindows -g lib/schMedA2.chromsize -w 100000 > lib/schMedA2.100kb.bed

rule cooltools_exp:
    input:
       rules.cooler_balance.output,
       cool = "04_cool/{GE}__{HIC}_exres_" + config["exres_cool"] + ".cool"
    output:
        touch("04_cool/done__cooltools_exp_{GE}__{HIC}"),
    log:
        "04_cool/log__cooltools_exp_{GE}__{HIC}"
    params: 
        opt = config["opt_cooltools_exp"],
        cis = "04_cool/{GE}__{HIC}_cis_exp.tsv",
        trans = "04_cool/{GE}__{HIC}_trans_exp.tsv",
        M = "lib/{GE}.1Mb_debug.bed",
        k = "lib/{GE}.100k.bed"
    threads: 1
    conda: 
        "../envs/cooler.yaml"
    shell:
        """
        cooltools compute-expected --regions {params.M} --contact-type cis -p {threads} {input.cool} -o {params.cis}
        cooltools compute-expected --regions {params.M} --contact-type trans -p {threads} {input.cool} -o {params.trans}
        """
