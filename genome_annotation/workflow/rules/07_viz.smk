###############################################################################

# dotPlotly

###############################################################################


rule self_mm2:
    input:
        assembly = ASSEMBLY
    output:
        touch("07_viz/done__01_paf"),
        all = "07_viz/01_paf/" + ASSEMBLY_ID + ".minimap2.paf",
        long = "07_viz/01_paf/" + ASSEMBLY_ID + "_long.minimap2.paf",
        long_assembly = "01_naming/" + ASSEMBLY_ID + "_larger_than_" + str(config["dotplot_minlen"]) + ".fasta"
    params:
        args = config["dotplotly_minimap2_args"],
        minlen = config["dotplot_minlen"]
    threads: 
        config["threads_minimap2"] 
    log:
        "logs/log__07_viz_01_self_mm2"
    conda:
        "../envs/dotplotly.yaml"
    shell:
        """
        minimap2 -t {threads} {params.args} \
        {input.assembly} {input.assembly} > {output.all} 2> {log}
        
        seqkit seq -m {params.minlen} {input.assembly} > {output.long_assembly} 2> {log}

        echo "short subset" >> {log}
        minimap2 -t {threads} {params.args} \
        {output.long_assembly} {output.long_assembly} > {output.long} 2>> {log}
        """


# TODO this has a messy mv command
rule dotplotly:
    input:
        all = "07_viz/01_paf/" + ASSEMBLY_ID + ".minimap2.paf",
        long = "07_viz/01_paf/" + ASSEMBLY_ID + "_long.minimap2.paf"
    output:
        touch("07_viz/done__02_dotplotly"),
        report("07_viz/02_dotplotly/" + ASSEMBLY_ID + ".minimap2.plot.png",
        caption="../report/07_02_dotplotly.rst", category="07 Vizualisation", subcategory="dotplotly"),
        report("07_viz/02_dotplotly/" + ASSEMBLY_ID + ".minimap2.plot.html",
        caption="../report/07_02_dotplotly.rst", category="07 Vizualisation", subcategory="dotplotly")
    params:
        args = config["dotplotly_args"],
        prefix_all =  ASSEMBLY_ID + ".minimap2.plot",
        prefix_long = ASSEMBLY_ID + "_long.minimap2.plot"
    threads:
        1 
    log:
        "logs/log__07_viz_02_dotplotly"
    conda:
        "../envs/dotplotly.yaml"
    shell:
        """
        Rscript workflow/scripts/pafCoordsDotPlotly.R -i {input.all} -o {params.prefix_all} {params.args} &> {log}
        Rscript workflow/scripts/pafCoordsDotPlotly.R -i {input.long} -o {params.prefix_long} {params.args} &>> {log}
        mv *plot* 07_viz/02_dotplotly/
        """


rule r07_viz_done:
    input:
        rules.dotplotly.output
    output:
        touch("07_viz/done__07_viz")
