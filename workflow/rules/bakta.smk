#%
# description: Annotate using Bakta - a comprehensive bacterial genome annotation tool.
# category: Functional Annotation
# link:
# - https://github.com/oschwengers/bakta
# references:
# - 'Schwengers O, et al. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. [Microb Genom. 2021 Nov;7(11):000685]'
#%

rule bakta_db_setup:
    output:
        db_file = "resources/GEMresources/bakta_db/db/rRNA.i1m",
        db_dir = directory("resources/GEMresources/bakta_db/db/amrfinderplus-db"),
    conda:
        "../envs/bakta.yaml"
    log: "logs/bakta/bakta_db_setup/bakta_db_setup.log"
    shell:
        """
        mkdir -p resources/GEMresources
        cd resources/GEMresources
        bakta_db download --output bakta_db --type full 2> {log}
        echo "Bakta database location: $(pwd)/bakta_db" >> {log}
        echo "Database download completed at: $(date)" >> {log}
        amrfinder_update --force_update --database /projects/panGEM/GEMresources/bakta_db/db/amrfinderplus-db
        """

rule extract_meta_bakta:
    input:
        fna = "data/interim/all/fasta/{strains_fna}.fna",
        samples_path = bgcflow_util_dir / "samples.csv",  
        assembly_report= "data/interim/all/assembly_report/{strains_fna}.json",
    output:
        org_info = "data/interim/all/bakta/{strains_fna}/organism_info.txt",
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/all/bakta/extract_meta_bakta/extract_meta_bakta-{strains_fna}.log"
    shell:
        """
        python workflow/bgcflow/bgcflow/data/get_organism_info.py {wildcards.strains_fna} \
            "{input.samples_path}" {input.assembly_report} data/interim/all/bakta/ 2>> {log}
        """

rule generate_locus_tags_bakta:
    """Generate unique locus tag prefix based on species name and accession"""
    input:
        org_info = "data/interim/all/bakta/{strains_fna}/organism_info.txt"
    output:
        locus_info = "data/interim/all/bakta/{strains_fna}/locus_info.json"
    params:
        strain = "{strains_fna}"
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/all/bakta/generate_locus_tags/generate_locus_tags-{strains_fna}.log"
    script:
        "../scripts/generate_locus_tags.py"
        
rule bakta:
    """Run Bakta annotation with computed locus tags and metadata"""
    input:
        fasta = "data/interim/all/fasta/{strains_fna}.fna",
        org_info = "data/interim/all/bakta/{strains_fna}/organism_info.txt",
        locus_info = "data/interim/all/bakta/{strains_fna}/locus_info.json",
        db_file = "resources/GEMresources/bakta_db/db/rRNA.i1m",
    output:
        gff = "data/interim/{stage}/bakta/{strains_fna}/{strains_fna}.gff3",
        faa = "data/interim/{stage}/bakta/{strains_fna}/{strains_fna}.faa",
        ffn = "data/interim/{stage}/bakta/{strains_fna}/{strains_fna}.ffn",
        gbff = "data/interim/{stage}/bakta/{strains_fna}/{strains_fna}.gbff",
        tsv = "data/interim/{stage}/bakta/{strains_fna}/{strains_fna}.tsv",
        txt = "data/interim/{stage}/bakta/{strains_fna}/{strains_fna}.txt",
    params:
        outdir = "data/interim/{stage}/bakta/{strains_fna}",
        prefix = "{strains_fna}",
        bakta_db = "resources/GEMresources/bakta_db/db",
    threads: 8
    conda:
        "../envs/bakta.yaml"
    log: "logs/{stage}/bakta/bakta/bakta-{strains_fna}.log"
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}

        # Extract genus, species from organism_info.txt (format: GENUS,SPECIES,STRAIN_ID)
        GENUS=$(cut -d "," -f 1 {input.org_info})
        SPECIES=$(cut -d "," -f 2 {input.org_info})

        # Extract locus tag from locus_info.json
        LOCUS_TAG=$(python3 -c "import json; data=json.load(open('{input.locus_info}')); print(data.get('locus_tag_prefix', ''))")

        # Log run parameters
        echo "=== Bakta Run Information ===" >> {log}
        echo "Timestamp: $(date)" >> {log}
        echo "Genome file: {input.fasta}" >> {log}
        echo "Genome ID: {wildcards.strains_fna}" >> {log}
        echo "Running Bakta for {input.fasta} with:" >> {log}
        echo "  --prefix '{params.prefix}'" >> {log}
        echo "  --locus-tag '$LOCUS_TAG'" >> {log}
        echo "  --genus '$GENUS'" >> {log}
        echo "  --species '$SPECIES'" >> {log}
        echo "  --threads {threads}" >> {log}
        echo "  --db {params.bakta_db}" >> {log}
        echo "===========================" >> {log}
        echo "" >> {log}

        # Run Bakta
        bakta \
            --db {params.bakta_db} \
            --force \
            --output {params.outdir} \
            --prefix {params.prefix} \
            --locus-tag "$LOCUS_TAG" \
            --genus "$GENUS" \
            --species "$SPECIES" \
            --threads {threads} \
            {input.fasta} &>> {log}

        echo "" >> {log}
        echo "=== Bakta Completed ===" >> {log}
        echo "Completion time: $(date)" >> {log}
        """

rule original_locus_tags_bakta:
    """Map original NCBI locus tags to bakta locus tags for all genomes in a species project"""
    input:
        ncbi_gff=fexpand("data/interim/all/gff/{accession}.gff", accession=RULE_FUNCTIONS.get("original_locus_tags", {}).get("accession", lambda w: [])),
        bakta_gbff=fexpand("data/interim/{{stage}}/bakta/{accession}/{accession}.gbff", accession=RULE_FUNCTIONS.get("original_locus_tags", {}).get("accession", lambda w: [])),
        samples="data/processed/species/samples/{name}.csv",
    output:
        csv="data/processed/{stage}/{name}/tables/df_locus_tag_mapping_bakta.csv",
    params:
        ncbi_gff="data/interim/all/gff/",
        bakta_gbff="data/interim/{stage}/bakta/",
    conda:
        "../envs/locus_tag_mapping.yaml"
    log:
        "logs/{stage}/original_locus_tags_bakta/original_locus_tags_bakta-{name}.log"
    shell:
        """
        python workflow/scripts/map_locus_tags_bakta.py \
            --samples {input.samples} \
            --ncbi_gff {params.ncbi_gff} \
            --bakta_gbff {params.bakta_gbff} \
            -o {output.csv} > {log} 2>&1
        """

rule compile_bakta_metadata:
    """Compile all per-genome bakta metadata into a single CSV file"""
    input:
        metadata_files = lambda wildcards: expand("data/interim/{{stage}}/bakta/{strains}/metadata.json",
                               strains=RULE_FUNCTIONS.get("bakta", {}).get(wildcards.stage, {}).get("samples", lambda: [])())
    output:
        csv = "data/processed/{stage}/{name}/tables/bakta_metadata.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log: "logs/{stage}/bakta/compile_bakta_metadata/compile_bakta_metadata-{name}.log"
    script:
        "../scripts/compile_bakta_metadata.py"
    