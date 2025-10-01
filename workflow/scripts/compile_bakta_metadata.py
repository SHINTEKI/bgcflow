#!/usr/bin/env python3
"""
Compile all per-genome bakta metadata into a single CSV file.
Adapted from bakta_pipeline.py for use in Snakemake workflow.
"""

import json
import csv
import os
from pathlib import Path

def compile_metadata_to_csv(metadata_files, output_csv):
    """
    Compile all metadata JSON files into a single CSV file.
    """
    metadata_list = []

    for metadata_file in metadata_files:
        if os.path.exists(metadata_file):
            try:
                with open(metadata_file, "r") as f:
                    data = json.load(f)
                    metadata_list.append(data)
                    print(f"Loaded metadata from {metadata_file}")
            except Exception as e:
                print(f"Error loading {metadata_file}: {e}")
        else:
            print(f"Warning: Metadata file not found: {metadata_file}")

    if metadata_list:
        # Get union of all keys from all metadata dictionaries
        fieldnames = sorted({k for md in metadata_list for k in md.keys()})

        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for md in metadata_list:
                writer.writerow(md)

        print(f"Compiled {len(metadata_list)} metadata entries into {output_csv}")
        return len(metadata_list)
    else:
        # Create empty CSV with basic headers
        with open(output_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["genome_id", "accession", "SpeciesName", "genus", "species"])

        print(f"No metadata found. Created empty CSV: {output_csv}")
        return 0

if __name__ == "__main__":
    # For Snakemake execution
    metadata_files = snakemake.input.metadata_files
    output_csv = snakemake.output.csv

    try:
        count = compile_metadata_to_csv(metadata_files, output_csv)
        print(f"Successfully compiled {count} metadata entries")

    except Exception as e:
        print(f"Error compiling metadata: {e}")
        raise