#%
# final_output: data/processed/{name}/roary/df_gene_presence_binary.csv
# description: Build pangenome using Roary.
# category: Comparative Genomics
# link:
# - https://github.com/sanger-pathogens/Roary
# references:
# - 'Andrew J. Page, Carla A. Cummins, Martin Hunt, Vanessa K. Wong, Sandra Reuter,
#   Matthew T. G. Holden, Maria Fookes, Daniel Falush, Jacqueline A. Keane, Julian
#   Parkhill, ''Roary: Rapid large-scale prokaryote pan genome analysis'', Bioinformatics,
#   2015;31(22):3691-3693 [doi:10.1093/bioinformatics/btv421](https://academic.oup.com/bioinformatics/article/31/22/3691/240757)'
#%
rule roary:
    input:
        gff=lambda wildcards: expand("data/interim/{stage}/prokka/{sample}/{sample}.gff", 
                                    stage=wildcards.stage, 
                                    sample=RULE_FUNCTIONS["roary"]["samples"](wildcards)),
        # gff=lambda wildcards: get_prokka_outputs(wildcards.name, filter_samples_qc(wildcards, get_samples_df())),
    output:
        core_alignment_header="data/interim/{stage}/roary/{name}/core_alignment_header.embl",
        pan_genome_reference="data/interim/{stage}/roary/{name}/pan_genome_reference.fa",
        accessory_header="data/interim/{stage}/roary/{name}/accessory.header.embl",
        clustered_proteins="data/interim/{stage}/roary/{name}/clustered_proteins",
        core_gene_alignment="data/interim/{stage}/roary/{name}/core_gene_alignment.aln",
        gene_presence_absence="data/interim/{stage}/roary/{name}/gene_presence_absence.csv",
        summary_statistics="data/interim/{stage}/roary/{name}/summary_statistics.txt",
        accessory_tab="data/interim/{stage}/roary/{name}/accessory.tab",
        accessory_binary_genes="data/interim/{stage}/roary/{name}/accessory_binary_genes.fa",
        core_accessory_header="data/interim/{stage}/roary/{name}/core_accessory.header.embl",
        numer_of_genes_in_pan_genome="data/interim/{stage}/roary/{name}/number_of_genes_in_pan_genome.Rtab",
        accessory_binary_genes_newick="data/interim/{stage}/roary/{name}/accessory_binary_genes.fa.newick",
        core_accessory="data/interim/{stage}/roary/{name}/core_accessory.tab",
        number_of_new_genes="data/interim/{stage}/roary/{name}/number_of_new_genes.Rtab",
        accessory_graph="data/interim/{stage}/roary/{name}/accessory_graph.dot",
        core_accessory_graph="data/interim/{stage}/roary/{name}/core_accessory_graph.dot",
        number_of_unique_genes="data/interim/{stage}/roary/{name}/number_of_unique_genes.Rtab",
    conda:
        "../envs/roary.yaml"
    params:
        i=80,
        g=lambda wildcards, output: max(80000, 40*len(RULE_FUNCTIONS["roary"]["samples"](wildcards))),
        roary_dir="data/interim/{stage}/roary/{name}/",
    threads: 64
    log:
        "logs/{stage}/roary/roary-{name}.log",
    shell:
        """
        rm -rf {params.roary_dir}
        roary -p {threads} -f {params.roary_dir} -i {params.i} -g {params.g} -e -n -r -v {input.gff} &>> {log}
        """

checkpoint roary_reassign_pangene_categories:
    input:
        gene_presence_binary="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        pangene_summary="data/interim/{stage}/roary/{name}/df_pangene_summary.csv"
    output:
        pangene="data/processed/{stage}/{name}/tables/df_roary_pangene_summary_reassigned.csv",
    log:
        "logs/{stage}/roary_reassign_pangene_categories/{name}.log"
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/scripts/alleleome_reassign_pangene_categories.py \
            --gp_binary {input.gene_presence_binary} \
            --summary {input.pangene_summary} \
            --output-file {output.pangene} 2>> {log}
        """

rule eggnog_roary:
    input:
        faa="data/interim/{stage}/roary/{name}/pan_genome_reference.fa",
        eggnog_db="resources/eggnog_db",
        dmnd="resources/eggnog_db/bacteria.dmnd",
    output:
        annotations="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations",
        xlsx="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations.xlsx",
        hits="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.hits",
        seed_orthologs="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.seed_orthologs",
        temp_dir=temp(directory("data/interim/{stage}/eggnog_roary/tmp/{name}"))
    params:
        eggnog_dir="data/interim/{stage}/eggnog_roary/{name}/",
    conda:
        "../envs/eggnog.yaml"
    threads: 8
    log:
        "logs/{stage}/eggnog-roary/eggnog-{name}.log",
    shell:
        """
        mkdir -p {output.temp_dir}
        emapper.py -i {input.faa} --translate --itype "CDS" --excel --cpu {threads} -o {wildcards.name} --output_dir {params.eggnog_dir} --data_dir {input.eggnog_db} --temp_dir {output.temp_dir} &>> {log}
        """

rule eggnog_roary_result_copy:
    input:
        annotations="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations",
        xlsx="data/interim/{stage}/eggnog_roary/{name}/{name}.emapper.annotations.xlsx",
    output:
        eggnog_xlsx="data/processed/{stage}/{name}/eggnog_roary/eggnog_roary.xlsx",
        eggnog_annotations="data/processed/{stage}/{name}/eggnog_roary/emapper.annotations"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/{stage}/eggnog-roary/eggnog-result-copy-{name}.log",
    shell:
        """
        cp {input.xlsx} {output.eggnog_xlsx}
        cp {input.annotations} {output.eggnog_annotations}
        """

rule deeptfactor_roary:
    input:
        pan_genome_reference="data/interim/{stage}/roary/{name}/pan_genome_reference.fa",
        resource="resources/deeptfactor/",
    output:
        df_deeptfactor="data/processed/{stage}/{name}/tables/df_deeptfactor_roary.tsv",
    conda:
        "../envs/deeptfactor.yaml"
    threads: 8
    log:
        "logs/{stage}/deeptfactor-roary/deeptfactor-roary-{name}.log",
    params:
        outdir="data/interim/{stage}/deeptfactor_roary/{name}",
    shell:
        """
        workdir=$PWD
        mkdir -p data/processed/{wildcards.name} 2>> {log}
        (cd {input.resource} && python tf_running.py \
            -i $workdir/{input.pan_genome_reference.fa} -o $workdir/{params.outdir} \
            -g cpu -cpu {threads}) 2>> {log}
        cp {params.outdir}/prediction_result.txt {output.df_deeptfactor} &>> {log}
        """

rule diamond_roary:
    input:
        faa="data/interim/{stage}/roary/{name}/pan_genome_reference.fa",
        resource="resources/deeptfactor/",
    output:
        diamond_interim="data/interim/{stage}/diamond/{name}/{name}_pangenome.dmnd",
        diamond_processed="data/processed/{stage}/{name}/diamond/{name}_pangenome.dmnd",
    conda:
        "../envs/antismash.yaml"
    threads: 8
    log:
        "logs/{stage}/diamond-roary/diamond-roary-{name}.log",
    shell:
        """
        diamond makedb --in {input.faa} -d {output.diamond_interim} -p {threads} &>> {log}
        cp {output.diamond_interim} {output.diamond_processed} &>> {log}
        """


rule roary_out:
    input:
        summary_statistics="data/interim/{stage}/roary/{name}/summary_statistics.txt",
        genomes_tree="data/processed/{stage}/{name}/automlst_wrapper/df_genomes_tree.csv",
        final_tree="data/processed/{stage}/{name}/automlst_wrapper/final.newick",
    output:
        gene_presence_binary="data/processed/{stage}/{name}/roary/df_gene_presence_binary.csv",
        gene_presence_locustag="data/processed/{stage}/{name}/roary/df_gene_presence_locustag.csv",
        summary="data/processed/{stage}/{name}/roary/df_pangene_summary.csv",
        summary_interim="data/interim/{stage}/roary/{name}/df_pangene_summary.csv",
        gene_presence_binary_interim="data/interim/{stage}/roary/{name}/df_gene_presence_binary.csv",
        gene_presence_locustag_interim="data/interim/{stage}/roary/{name}/df_gene_presence_locustag.csv",
    params:
        roary_interim_dir="data/interim/{stage}/roary/{name}/",
        roary_processed_dir="data/processed/{stage}/{name}/roary",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/{stage}/roary/roary-out-{name}.log",
    shell:
        """
        python workflow/bgcflow/bgcflow/data/make_pangenome_dataset.py {params.roary_interim_dir} {params.roary_processed_dir} {input.final_tree} {input.genomes_tree} 2>> {log}
        """
