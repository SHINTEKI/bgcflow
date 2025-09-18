rule custom_samples:
    output:
        samples_file="data/interim/all/custom_genomes/samples.csv",
    run:
        if CUSTOM_SAMPLES is None:
            dummy_df = pd.DataFrame(data=[], columns=["source", "path"])
            dummy_df.index.name = "genome_id"
            dummy_df.to_csv(output.samples_file)
        else:
            CUSTOM_SAMPLES[["source", "path"]].to_csv(output.samples_file)