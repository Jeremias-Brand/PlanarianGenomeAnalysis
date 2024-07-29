###############################################################################

# dotPlotly

###############################################################################

rule comp_mm2:
    input:
        GE1 = config['genomeDIR'] + "{GE1}" + config['genomeEXT'],
        GE2 = config['genomeDIR'] + "{GE2}" + config['genomeEXT']
    output:
        touch("dotplots/done__{GE1}__{GE2}_paf"),
        all = "dotplots/01_paf/{GE1}__{GE2}.minimap2.paf"
        #long = "dotplots/01_paf/{GE1}__{GE2}_long.minimap2.paf"
    params:
        args = config["dotplotly_minimap2_args"],
        minlen = config["dotplot_minlen"],
        long_GE1 = "01_naming/{GE1}_larger_than_" + str(config["dotplot_minlen"]) + ".fa",
        long_GE2 = "01_naming/{GE2}_larger_than_" + str(config["dotplot_minlen"]) + ".fa"
    threads: 
        config["threads_minimap2"] 
    log:
        "logs/log__dotplots_comp_mm2_{GE1}__{GE2}"
    conda:
        "../envs/dotplotly.yaml"
    shell:
        """
        mkdir -p 01_naming
        #echo "generate short subset" > {log}
        #seqkit seq -m {params.minlen} {input.GE1} > {params.long_GE1} 2>> {log}
        #seqkit seq -m {params.minlen} {input.GE2} > {params.long_GE2} 2>> {log}

        echo "Starting minimap2 alignment 1" >> {log}
        minimap2 -t {threads} {params.args} \
        {input.GE1} {input.GE2} > {output.all} 2>> {log}
        """


#rule self_mm2:
#    input:
#        assembly = ASSEMBLY
#    output:
#        touch("07_viz/done__01_paf"),
#        all = "07_viz/01_paf/" + ASSEMBLY_ID + ".minimap2.paf",
#        long = "07_viz/01_paf/" + ASSEMBLY_ID + "_long.minimap2.paf",
#        long_assembly = "01_naming/" + ASSEMBLY_ID + "_larger_than_" + str(config["dotplot_minlen"]) + ".fasta"
#    params:
#        args = config["dotplotly_minimap2_args"],
#        minlen = config["dotplot_minlen"]
#    threads: 
#        config["threads_minimap2"] 
#    log:
#        "logs/log__07_viz_01_self_mm2"
#    conda:
#        "../envs/dotplotly.yaml"
#    shell:
#        """
#        minimap2 -t {threads} {params.args} \
#        {input.assembly} {input.assembly} > {output.all} 2> {log}
#        
#        seqkit seq -m {params.minlen} {input.assembly} > {output.long_assembly} 2> {log}
#
#        echo "short subset" >> {log}
#        minimap2 -t {threads} {params.args} \
#        {output.long_assembly} {output.long_assembly} > {output.long} 2>> {log}
#        """
#
#
## TODO this has a messy mv command
rule dotplotly:
    input:
        all = "dotplots/01_paf/{GE1}__{GE2}.minimap2.paf"
        # long = "dotplots/01_paf/{GE1}__{GE2}_long.minimap2.paf"
    output:
        touch("dotplots/done__02_dotplotly_{GE1}__{GE2}"),
        "dotplots/02_dotplotly/{GE1}__{GE2}.minimap2.plot.png",
        #caption="../report/07_02_dotplotly.rst", category="07 Vizualisation", subcategory="dotplotly"),
        "dotplots/02_dotplotly/{GE1}__{GE2}.minimap2.plot.html"
        #caption="../report/07_02_dotplotly.rst", category="07 Vizualisation", subcategory="dotplotly")
    params:
        args = config["dotplotly_args"],
        prefix_all =  ".minimap2.plot",
        prefix_long = "_long.minimap2.plot"
    threads:
        1 
    log:
        "logs/log__dotplots_02_dotplotly_{GE1}__{GE2}"
    conda:
        "../envs/dotplotly.yaml"
    shell:
        """
        Rscript workflow/scripts/pafCoordsDotPlotly.R -i {input.all} -o {wildcards.GE1}__{wildcards.GE2}{params.prefix_all} {params.args} &> {log}
        mv *plot*.html dotplots/02_dotplotly/
        mv *plot*.png dotplots/02_dotplotly/
        """
#
#
#
