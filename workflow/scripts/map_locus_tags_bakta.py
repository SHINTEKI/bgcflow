#!/usr/bin/env python3
"""
Map locus tags for Bakta annotations.
Wrapper around map_locus_tags.py that handles GBFF files.
"""

import os
import json
from pathlib import Path
import argparse
from BCBio import GFF
from Bio import SeqIO
import pandas as pd
from sortedcontainers import SortedDict

def initialize_parser(parser):
    parser.description = "Map original NCBI locus tags to Bakta locus tags."
    parser.add_argument(
        "--samples",
        type=str,
        help="Table with genome accessions under the column 'genome_id'.",
    )
    parser.add_argument(
        "--ncbi_gff",
        type=str,
        required=False,
        default=None,
        help="GFF file from NCBI (with the original locus_tags).",
    )
    parser.add_argument(
        "--bakta_gbff",
        type=str,
        required=True,
        help="Bakta GBFF file directory (with the bakta locus_tags).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file with mapping.",
    )

def locus_tag_mapping_bakta(
    samples_path,
    ncbi_gff_path,
    bakta_gbff_path,
    output_path,
):
    try:
        accessions = pd.read_csv(samples_path, usecols=["genome_id"])
        accessions = accessions["genome_id"].to_list()
    except pd.errors.EmptyDataError:
        accessions = []

    bakta_gbff_path = Path(bakta_gbff_path)
    ncbi_gff_path = Path(ncbi_gff_path) if ncbi_gff_path else None

    df = None
    with open(output_path, "w") as f_output:
        for i_accession, accession in enumerate(accessions):
            print(f"Accession {i_accession + 1}/{len(accessions)} ({accession})", flush=True)
            mapping = []
            feature_dict = SortedDict()

            # Parse NCBI GFF for original locus tags
            if ncbi_gff_path:
                accession_file = ncbi_gff_path / f"{accession}.gff"
                if os.path.exists(accession_file) and os.stat(accession_file).st_size > 0:
                    with open(accession_file, "r") as f:
                        for record in GFF.parse(f):
                            for feature in record.features:
                                try:
                                    feature_type = feature.type
                                    if not feature_type in ["gene", "pseudogene"]:
                                        continue
                                    start_pos = int(feature.location.start)
                                    end_pos = int(feature.location.end)
                                    strand = feature.location.strand
                                    locus_tag = feature.qualifiers["locus_tag"][0]
                                    gene_key = "gene" if "gene" in feature.qualifiers else "Name"
                                    gene = feature.qualifiers[gene_key][0]
                                    synonyms = set()
                                    if len(feature.qualifiers[gene_key]) > 1:
                                        for synonym in feature.qualifiers[gene_key][1:]:
                                            synonyms.add(synonym)
                                    for synonym in feature.qualifiers.get("gene_synonym", []):
                                        synonyms.add(synonym)
                                    synonyms = list(synonyms)
                                    feature_key = (start_pos, end_pos, strand)

                                    feature_present = False
                                    if feature_key[0] in feature_dict:
                                        for f in feature_dict[feature_key[0]]:
                                            if f[0] == feature_key:
                                                feature_present = True
                                                break
                                    if not feature_present:
                                        if not feature_key[0] in feature_dict:
                                            feature_dict[feature_key[0]] = []
                                        feature_dict[feature_key[0]].append((feature_key, (locus_tag, feature_type, gene, synonyms)))
                                except Exception as e:
                                    print(f"Could not properly parse feature {feature.id}: {e}")
                else:
                    print(f"Skipping, because file {accession_file} is empty or doesn't exist.")

            print("Built feature_dict", flush=True)

            # Parse Bakta GBFF file - try both .gbff and .gbk extensions
            bakta_file = bakta_gbff_path / f"{accession}.gbff"
            if not bakta_file.exists():
                bakta_file = bakta_gbff_path / f"{accession}.gbk"

            if bakta_file.exists():
                for record in SeqIO.parse(bakta_file, "genbank"):
                    for i, feature in enumerate(record.features):
                        feature_type = feature.type
                        if not feature_type == "CDS":
                            continue
                        start_pos = int(feature.location.start)
                        end_pos = int(feature.location.end)
                        strand = feature.location.strand
                        feature_key = (start_pos, end_pos, strand)
                        locus_tag = feature.qualifiers["locus_tag"][0]
                        gene = feature.qualifiers.get("gene", [None])[0]

                        original_info = None
                        exact_match = False
                        for f in feature_dict.get(feature_key[0], []):
                            if f[0] == feature_key:
                                original_info = f[1]
                                exact_match = True
                                break
                        if not exact_match:
                            min_diff = None
                            for sp in feature_dict.irange(start_pos - 300, start_pos + 300):
                                for f in feature_dict[sp]:
                                    k, v = f
                                    if k[2] != strand:
                                        continue
                                    f_start_pos = k[0]
                                    f_end_pos = k[1]
                                    start_diff = abs(f_start_pos - start_pos)
                                    end_diff = abs(f_end_pos - end_pos)
                                    if (end_pos - f_start_pos) > 30 and (f_end_pos - f_start_pos) > 30 and (start_diff % 3) == 0 and (end_diff % 3) == 0 and ((start_diff < 210 and end_diff < 210) or (min(start_diff, end_diff) < 21 and max(start_diff, end_diff) < 300)):
                                        if min_diff is None or min_diff >= (start_diff + end_diff):
                                            feature_key = k
                                            original_info = v

                        if original_info is not None:
                            mapping.append(
                                {
                                    "bakta_locus_tag": locus_tag,
                                    "genome": accession,
                                    "bakta_gene": gene,
                                    "bakta_start_pos": start_pos,
                                    "bakta_end_pos": end_pos,
                                    "bakta_strand": strand,
                                    "original_locus_tag": original_info[0],
                                    "original_type": original_info[1],
                                    "original_gene": original_info[2],
                                    "original_synonyms": ",".join(original_info[3]),
                                    "original_start_pos": feature_key[0],
                                    "original_end_pos": feature_key[1],
                                    "original_strand": feature_key[2],
                                    "exact_match": exact_match,
                                }
                            )
            else:
                print(f"Warning: Bakta file not found for {accession}")

            first_write = (df is None)

            df = pd.DataFrame.from_records(mapping)
            if len(df) > 0:
                if first_write:
                    df.to_csv(f_output, index=False)
                else:
                    df.to_csv(f_output, index=False, header=False)

def run(args):
    locus_tag_mapping_bakta(
        args.samples,
        args.ncbi_gff,
        args.bakta_gbff,
        args.output,
    )

def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()