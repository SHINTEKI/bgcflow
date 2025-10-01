#!/usr/bin/env python3
"""
Generate locus tag prefix based on species name and accession.
Reads from organism_info.txt (format: GENUS,SPECIES,STRAIN_ID)
"""

import json
import re
import os

def generate_locus_tag_prefix(genus, species, strain_name):
    """
    Generate unique locus tag prefix based on genus, species and accession.
    """
    # Create species prefix from genus + species (first 3 letters each)
    species_name = f"{genus} {species}"
    words = species_name.split()
    species_prefix = "".join([word[:3] for word in words]).lower()

    # Extract accession from strain name
    match = re.match(r"(GC[AF]_\d+\.\d+)", strain_name)
    if not match:
        raise ValueError("Accession not found in strain name for locus tag generation.")

    accession = match.group(1)
    part = accession.split('.')[0]

    prefix = (species_prefix + part).upper()
    print(f"Generated locus tag prefix: {prefix}")

    return prefix

if __name__ == "__main__":
    # For Snakemake execution
    input_org_info = snakemake.input.org_info
    output_locus_info = snakemake.output.locus_info
    strain_name = snakemake.params.strain

    try:
        # Load organism info from organism_info.txt
        # Format: GENUS,SPECIES,STRAIN_ID
        with open(input_org_info, "r") as f:
            org_data = f.read().strip().split(",")

        genus = org_data[0] if len(org_data) > 0 else ""
        species = org_data[1] if len(org_data) > 1 else ""
        strain_id = org_data[2] if len(org_data) > 2 else ""

        # Generate locus tag prefix
        locus_tag_prefix = generate_locus_tag_prefix(genus, species, strain_name)

        # Create locus info dictionary
        locus_info = {
            "strain": strain_name,
            "locus_tag_prefix": locus_tag_prefix,
            "genus": genus,
            "species": f"{genus} {species}",
            "strain_id": strain_id
        }

        # Save locus info to JSON file
        with open(output_locus_info, "w") as f:
            json.dump(locus_info, f, indent=4)

        print(f"Locus tag info saved to {output_locus_info}")

    except Exception as e:
        print(f"Error generating locus tag for {strain_name}: {e}")
        raise