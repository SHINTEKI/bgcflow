checkpoint qc:
    input:
        seqfu=fexpand("data/processed/{stage}/{name}/tables/df_seqfu_stats.csv", name=RULE_FUNCTIONS["qc"]["names"], stage=RULE_FUNCTIONS["qc"]["stages"]),
        checkm=fexpand("data/processed/{stage}/{name}/tables/df_combined_stats.csv", name=RULE_FUNCTIONS["qc"]["names"], stage=RULE_FUNCTIONS["qc"]["stages"]),
    output: "data/processed/qc/qc_passed.csv",
    run:
        import pandas as pd
        df_seqfu = None
        for f in input.seqfu:
            df_s = pd.read_csv(f, sep=",", index_col=0, header=0)
            df_seqfu = pd.concat((df_seqfu, df_s))
        df_checkm = None
        for f in input.checkm:
            df_s = pd.read_csv(f, sep=",", index_col=0, header=0)
            df_checkm = pd.concat((df_checkm, df_s))
        df = df_seqfu.merge(df_checkm, how='outer', left_index=True, right_index=True)
        df["passed"] = ((df["N50"] > 50_000) & (df["# contigs"] < 200) & (df["Completeness"] > 95) & (df["Contamination"] < 5))
        df.to_csv(output[0], columns=["passed"], header=True, index=True, index_label="genome_id")
