import requests
import multiprocessing
import time
from pathlib import Path

# Read release version from config
try:
    gtdb_release = config["rule_parameters"]["install_gtdbtk"]["release"]
    gtdb_release_version = config["rule_parameters"]["install_gtdbtk"][
        "release_version"
    ]
except KeyError:
    gtdb_release = "220.0"
    gtdb_release_version = "220.0"
sys.stderr.write(f"Checking GTDB API...\n")
sys.stderr.write(f" - GTDB API | Grabbing metadata using GTDB release version: {gtdb_release_version}\n")

# Testing connection to GTDB API
gtdb_api_url = "https://gtdb-api.ecogenomic.org/status/db"

gtdb_offline_mode = False  # Flag to indicate if the GTDB API is offline

if "rule_parameters" in config.keys():
    if "use_gtdb_api" in config["rule_parameters"]:
        if config["rule_parameters"]["use_gtdb_api"] == False:
            gtdb_offline_mode = True

if not gtdb_offline_mode:
    try:
        sys.stderr.write(f" - GTDB API | Testing connection to: {gtdb_api_url}\n")
        response = requests.get(gtdb_api_url)
        response.raise_for_status()  # Raise an exception for HTTP errors (4xx, 5xx)
        data = response.json()  # Assuming the API returns JSON data
        sys.stderr.write(f' - GTDB API | Database is online: {data["online"]}\n')

    except requests.exceptions.RequestException as e:
        sys.stderr.write(f" - GTDB API | Error: {e}\n")
        sys.stderr.write(" - GTDB API | Failed to connect to the GTDB API.\n")
        sys.stderr.write(" - GTDB API | It is possible to continue in offline mode. This will return empty taxonomic information for all NCBI genomes!\n")

        # Check if the error is due to a 504 or 404 status code
        if response.status_code == 504 or response.status_code == 404:
            countdown_seconds = 30
            while countdown_seconds > 0:
                sys.stderr.write(f"\rDo you want to continue in offline mode? (yes/no/stop) (Time left: {countdown_seconds} seconds): ")
                sys.stderr.flush()
                time.sleep(1)  # Wait for 1 second
                countdown_seconds -= 1

            sys.stderr.write("\rDo you want to continue in offline mode? (yes/no/stop) (Time left: 0 seconds): \n")
            sys.stderr.flush()

            user_input = get_user_input_with_timeout("", 0)  # Get user input with no timeout
            if user_input is not None and user_input.strip().lower() == 'yes':
                gtdb_offline_mode = True  # Continue in offline mode
            elif user_input is not None and user_input.strip().lower() == 'stop':
                raise
            else:
                # timed out, continue anyway
                #sys.stderr.write("WARNING: No response, continuing BGCFlow in online mode anyway...")
                pass

else:
    sys.stderr.write("WARNING! GTDB API offline mode enabled. This will return empty taxonomic information for all NCBI genomes!\n")

sys.stderr.write(f" - GTDB API | Searching in offline mode: {gtdb_offline_mode}\n")

rule gtdb_install_table:
    output:
        table=f"resources/gtdb_download/bac120_metadata_r{str(gtdb_release).split('.')[0]}.tsv"
    priority: 50
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/gtdb/gtdb_install_table/gtdb_install_table.log",
    params:
        gtdb_link_au=f"https://data.ace.uq.edu.au/public/gtdb/data/releases/release{str(gtdb_release).split('.')[0]}/{gtdb_release}/bac120_metadata_r{str(gtdb_release).split('.')[0]}.tsv.gz",
        gtdb_link_main=f"https://data.gtdb.ecogenomic.org/releases/release{str(gtdb_release).split('.')[0]}/{gtdb_release}/bac120_metadata_r{str(gtdb_release).split('.')[0]}.tsv.gz",
        table_gz=lambda wildcards, output: f"{output.table}.gz",
        download_dir=lambda wildcards, output: Path(output.table).parent,
    shell:
        """
            # Check if the gz file already exists and is valid
            if [ -f "{params.table_gz}" ] && gzip -t "{params.table_gz}" 2>/dev/null; then
                echo "Found valid GTDB table file, skipping download" >> {log}
            else
                echo "Downloading GTDB table..." >> {log}
                # Try Australian server first (faster), fallback to main server
                echo "Attempting download from Australian server (data.ace.uq.edu.au)..." >> {log}
                if aria2c -x 16 -s 16 -c -d {params.download_dir} {params.gtdb_link_au} &>> {log}; then
                    echo "Successfully downloaded from Australian server" >> {log}
                else
                    echo "Australian server failed, trying main server (data.gtdb.ecogenomic.org)..." >> {log}
                    aria2c -x 16 -s 16 -c -d {params.download_dir} {params.gtdb_link_main} &>> {log} && echo "Successfully downloaded from main server" >> {log}
                fi
            fi
            # skip gunzip to save time
            gunzip -c '{params.table_gz}' > {output.table} 2>> {log}
            # Clean up aria2 control file if exists, but keep the tar.gz
            rm -f '{params.table_gz}.aria2' 2>/dev/null || true
            echo "GTDB table setup complete. Keeping gz file for future use." >> {log}
        """

rule gtdb_prep:
    input:
        table=f"resources/gtdb_download/bac120_metadata_r{str(gtdb_release).split('.')[0]}.tsv",
        samples_csv="data/interim/{stage}/ncbi_datasets/{taxon}.csv",
    output:
        gtdb_jsonl="data/interim/{stage}/gtdb/{taxon}.jsonl",
    log:
        "logs/{stage}/gtdb/gtdb_prep/gtdb_prep-{taxon}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    params:
        # samples_path=bgcflow_util_dir / "samples.csv",
        gtdb_paths=GTDB_PATHS,
        version=f"R{str(gtdb_release).split('.')[0]}",
        gtdb_table = f"bac120_metadata_r{str(gtdb_release).split('.')[0]}.tsv",
        offline=gtdb_offline_mode,
        api_base="https://gtdb-api.ecogenomic.org",
    shell:
        """
            python workflow/bgcflow/bgcflow/data/gtdb_prep_from_table_or_api.py {input.samples_csv} {input.table} {params.version} {output.gtdb_jsonl} '{params.gtdb_paths}' {params.offline} {params.api_base} > {log} 2>&1
        """

rule gtdb_jsonl_combine:
    input:
        gtdb_jsonl=lambda wildcards: expand("data/interim/{{stage}}/gtdb/{taxon}.jsonl", taxon=RULE_FUNCTIONS["gtdb_improved"]["taxons"]()),
    output:
        gtdb_jsonl="data/interim/{stage}/gtdb/all/gtdb.jsonl",
    shell:
        """
            touch {output.gtdb_jsonl}
            for f in {input.gtdb_jsonl}
            do
                cat $f >> {output.gtdb_jsonl}
            done
        """
rule gtdb_jsonl_to_json:
    input:
        gtdb_jsonl="data/interim/{stage}/gtdb/all/gtdb.jsonl",
    output:
        gtdb_json="data/interim/{stage}/gtdb/{accession}.json",
    shell:
        """
            grep '^{{"genome_id": "{wildcards.accession}"' {input.gtdb_jsonl} > {output.gtdb_json}
        """


checkpoint fix_gtdb_taxonomy:
    input:
        gtdb_jsonl="data/interim/{stage}/gtdb/{taxon}.jsonl",
    output:
        meta="data/interim/{stage}/gtdb/{taxon}/tables/df_gtdb_meta.csv",
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/{stage}/gtdb/fix_gtdb_taxonomy/fix_gtdb_taxonomy-{taxon}.log",
    priority: 50
    shell:
        """
            python workflow/bgcflow/bgcflow/data/fix_gtdb_taxonomy.py {input.gtdb_jsonl} {output.meta} 2> {log}
        """
