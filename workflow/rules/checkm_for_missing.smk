CHECKM_COLMAPPING = {
    "CheckM completeness": "Completeness",
    "CheckM contamination": "Contamination",
    "Assembly Stats GC Percent": "GC",
    # "Assembly Stats Number of Scaffolds": "# scaffolds",
    "Assembly Stats Number of Contigs": "# contigs",
    # "Assembly Stats Total Sequence Length": "Genome size",
    # "Assembly Stats Contig N50": "N50 (contigs)",
    # "Assembly Stats Scaffold N50": "N50 (scaffolds)",
}

if config.get("use_ncbi_data_for_checkm", False):
    checkpoint checkm_find_missing_data:
        input:
            overview_csv="data/interim/{stage}/ncbi_datasets/{taxon}.csv",
        output:
            missing="data/interim/{stage}/checkm_missing/{taxon}/{taxon}_missing.txt",
        log:
            "logs/{stage}/checkm/checkm_missing_{taxon}.log",
        run:
            import pandas as pd
            df = pd.read_csv(input.overview_csv, low_memory=False, index_col=0, header=0)
            if not all((col in df.columns) for col in CHECKM_COLMAPPING.keys()):
                missing = df.index.to_series()
            else:
                missing = df.index[pd.isna(df[list(CHECKM_COLMAPPING.keys())]).any(axis="columns")].to_series()
            missing.to_csv(output.missing, index=False, header=False)

    def get_checkm_missing_fasta_inputs_for_name(stage, name):
        try:
            accessions = pd.read_csv(checkpoints.checkm_find_missing_data.get(stage=stage, taxon=name).output.missing, header=None, index_col=False)
        except pd.errors.EmptyDataError:
            return []
        if accessions.empty:
            return []
        accessions = accessions[accessions.columns[0]].to_list()
        return [f"data/interim/all/fasta/{s}.fna" for s in accessions]

    rule checkm_for_missing:
        input:
            fna=lambda wildcards: get_checkm_missing_fasta_inputs_for_name(wildcards.stage, wildcards.name),
            checkm_db="resources/checkm/",
            missing="data/interim/{stage}/checkm_missing/{name}/{name}_missing.txt",
        output:
            tempdir=temp(directory("data/interim/{stage}/checkm_missing_temp/{name}/{name}_fna")),
            fna=temp(directory("data/interim/{stage}/checkm_missing/{name}/{name}_fna")),
            stat="data/interim/{stage}/checkm_missing/{name}/{name}/storage/bin_stats_ext.tsv",
            checkm_dir=directory("data/interim/{stage}/checkm_missing/{name}/{name}"),
        conda:
            "../envs/checkm.yaml"
        log:
            "logs/{stage}/checkm/checkm_{name}.log",
        params:
            checkm_log="data/interim/{stage}/checkm/{name}/checkm_missing_{name}.log",
        threads: 16
        shell:
            """
            if [ -s {input.missing} ]; then
                mkdir -p {output.fna}
                mkdir -p {output.tempdir}
                mkdir -p {output.checkm_dir}/storage
                mkdir -p {output.checkm_dir}/bins
                for f in {input.fna}; do cp $f {output.fna}/.; done
                TEMPDIR={output.tempdir} checkm lineage_wf -t {threads} --reduced_tree -x fna {output.fna} {output.checkm_dir} &>> {log}
            else
                mkdir -p {output.checkm_dir}
                mkdir -p {output.checkm_dir}/storage
                mkdir -p {output.checkm_dir}/bins
                mkdir -p {output.tempdir}
                touch {output.stat}
                mkdir -p {output.fna}
            fi
            """

    rule checkm_missing_out:
        input:
            stat="data/interim/{stage}/checkm_missing/{name}/{name}/storage/bin_stats_ext.tsv",
        output:
            stat_processed="data/interim/{stage}/checkm_missing/{name}/{name}_stats.csv",
        log:
            "logs/{stage}/checkm/checkm_missing_out_{name}.log",
        params:
            checkm_json=directory("data/interim/{stage}/checkm_missing/{name}/{name}/json/"),
        conda:
            "../envs/bgc_analytics.yaml"
        shell:
            """
            if [ -s {input.stat} ]; then
                mkdir -p {params.checkm_json}
                python workflow/bgcflow/bgcflow/data/get_checkm_data.py {input.stat} {params.checkm_json} {output.stat_processed} 2>> {log}
            else
                touch {output.stat_processed}
            fi
            """

    rule checkm_missing_combine:
        input:
            checkm_stats="data/interim/{stage}/checkm_missing/{name}/{name}_stats.csv",
            overview_csv="data/interim/{stage}/ncbi_datasets/{name}.csv",
        output:
            stats="data/processed/{stage}/{name}/tables/df_combined_stats.csv",
        log:
            "logs/{stage}/checkm/checkm_missing_combine_{name}.log",
        run:
            import pandas as pd
            df = pd.read_csv(input.overview_csv, low_memory=False, index_col=0, header=0)
            try:
                df_checkm = pd.read_csv(input.checkm_stats, low_memory=False, index_col=0, header=0)
            except:
                df_checkm = None
            if all((col in df.columns) for col in CHECKM_COLMAPPING.keys()):
                df_combined = df[list(CHECKM_COLMAPPING.keys())].copy()
                df_combined.rename(columns=CHECKM_COLMAPPING, inplace=True)

                # Correct for percentages vs. fractions
                df_combined["GC"] = df_combined["GC"]/100.0

                if not df_checkm is None:
                    df_combined.update(df_checkm)
            else:
                df_combined = df_checkm[list(CHECKM_COLMAPPING.values())]
            df_combined.to_csv(output.stats)

else:
    rule use_checkm:
        input:
            checkm_stats="data/processed/{stage}/{name}/tables/df_checkm_stats.csv",
        output:
            stats="data/processed/{stage}/{name}/tables/df_combined_stats.csv",
        log:
            "logs/{stage}/checkm/checkm_missing_use_checkm_{name}.log",
        run:
            import pandas as pd
            df_checkm[list(CHECKM_COLMAPPING.values())].to_csv(output.stats)
