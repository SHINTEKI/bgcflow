import json
import logging
import sys
from pathlib import Path

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def assess_gtdb_json_file(jsonl):
    """
    Given a jsonl file, assess whether the each accession can be found via GTDB-API or not.
    Return a genome_id list when accession's taxonomy cannot be found.
    """
    logging.info(f"Assessing genome lists")

    genomeid_list = []
    with open(jsonl, "r", encoding="utf-8") as f:
        for line in f:
            record = json.loads(line)
            genome_id = record.get("genome_id")
            # assess if the required keys are present
            try:
                gtdb_release = record["gtdb_release"]
                metadata = record["metadata"]
                if "Genome not found" in metadata["detail"]:
                    logging.debug(
                    f" - {genome_id} : {metadata['detail']} in GTDB-API release {gtdb_release}"
                    )
                    genomeid_list.append(genome_id)
                elif type(metadata["genome"]["accession"]) == str:
                    logging.debug(
                    f" - {genome_id} can be found via GTDB-API release {gtdb_release}"
                    )
                    pass

            except KeyError:
                logging.debug(f" - {genome_id} does not have metadata")
                genomeid_list.append(genome_id)
    return genomeid_list
    # with open(item, "r") as json_file:
        # content = json_file.read().strip()
        
        # # Handle empty or whitespace-only files
        # if not content:
        #     genome_id = Path(item).stem
        #     logging.debug(f" - {genome_id} : Empty JSON file, needs GTDB-Tk classification")
        #     return genome_id
        # # Handle non-empty files 
        # data = json.loads(content)       
        # genome_id = data["currentAccession"]
        # try:
        #     gtdb_release = data["gtdb_release"]
        #     metadata = data["metadata"]
        #     if "Genome not found" in metadata["detail"]:
        #         logging.debug(
        #             f" - {genome_id} : {metadata['detail']} in GTDB-API release {gtdb_release}"
        #         )
        #         return genome_id
        #     elif type(metadata["genome"]["accession"]) == str:
        #         logging.debug(
        #             f" - {genome_id} can be found via GTDB-API release {gtdb_release}"
        #         )
        #         return None

        # except KeyError:
        #     logging.debug(f" - {genome_id} does not have metadata")
        #     return genome_id


def generate_symlink_gtdbtk(fna_paths_list, genome_id_list, outdir):
    """
    Given a fna paths list and a genome id list, generate a symlink to a desired location
    """

    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    index = {Path(p).stem: Path(p).resolve() for p in fna_paths_list}
    
    for gid in genome_id_list:
        gid = str(gid)
        if gid not in index:
            logging.debug(f"[ERROR] Matched fna file not found: {gid}", file=sys.stderr)
            sys.exit(1)

        src = index[gid]
        dst = outdir / src.name 
        dst.symlink_to(src) 
        logging.info(f"[OK] {dst} -> {src}")

    # if genome_id is not None:
    #     outfile = outdir / f"{genome_id}.fna"
    #     logging.info(f"Generating input files for GTDB-tk: {outfile}")
    #     outfile.symlink_to(input_fna)
        # logging.info(f"Symlink from {input_fna} to {outfile} generated")
        # return outfile.stem  # Return the actual filename stem used
    # else:
    #     return None


def input_handling(input_list, category, suffix=".json"):
    input_list = Path(input_list)
    if input_list.is_file() and input_list.suffix == suffix:
        logging.info(f"Getting {category} from a single file: {input_list}")
        input_list_files = input_list

    elif input_list.is_file() and input_list.suffix == ".txt":
        logging.info(f"Getting {category} from a text file: {input_list}")
        with open(input_list, "r") as file:
            file_content = [i.strip("\n") for i in file.readlines()]
            if len(file_content) == 1:
                # Paths space-separated on a single line
                logging.info(
                    " - Detecting space-separated input in a single line format."
                )
                paths = file_content[0].split()
            else:
                # Paths written on separate lines
                logging.info(" - Detecting input in a multi-line format.")
                paths = file_content
            input_list_files = [
                Path(path) for path in paths if Path(path).suffix == suffix
            ]
        logging.info(f"Getting {len(paths)} paths to handle")
    else:
        input_list_files = [
            Path(file)
            for file in str(input_list).split()
            if Path(file).suffix == suffix
        ]
        logging.info(
            f"Getting {category} from the given list of {len(input_list_files)} files..."
        )
    return input_list_files


def gtdbtk_prep(fna_list, gtdb_jsonl, outdir, output_txt):
    """
    Given a list fna paths and a gtdb jsonl, generate a symlinks to a desired location
    if genome_id cannot be found via GTDB API or gtdbtk table
    """
    fna_paths_list = input_handling(fna_list, "fna files", suffix=".fna")
    # Ensure fna files and gtdb_jsonl have the same length
    with open(gtdb_jsonl, "r", encoding="utf-8") as f:
        count = sum(1 for _ in f)
    logging.info(f"Number of fna files: {len(fna_paths_list)}")
    logging.info(f"Number of gtdb json lines: {count}")
    assert len(fna_paths_list) == count, "Number of fna files and gbtk jsons don't equal"

    # assess whether genome_id has been classified via GTDB-API or table, return a list of genome_id to be processeed
    genome_id_list = assess_gtdb_json_file(gtdb_jsonl)

    generate_symlink_gtdbtk(fna_paths_list, genome_id_list, str(outdir))
    
        # if fnafile is not None:
        #     logging.debug(f"Adding {fnafile} to input list for gtdbtk")
        #     input_list.append(fnafile)
    
    
    # generate_symlink_gtdbtk(str(input_fna[0]), str(gtdb_json), str(outdir))
    #     logging.debug(f"fnafile: {fnafile}")
    #     if fnafile is not None:
    #         logging.debug(f"Adding {fnafile} to input list for gtdbtk")
    #         input_list.append(fnafile)
    with open(output_txt, "w") as f:
        for item in genome_id_list:
            f.write("%s\n" % item)
    logging.info(f"Successfully processed {len(genome_id_list)} genomes for GTDB-Tk")
    return


if __name__ == "__main__":
    try:
        gtdbtk_prep(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
        logging.info("gtdbtk_prep completed successfully")
        sys.exit(0)
    except Exception as e:
        logging.error(f"gtdbtk_prep failed with error: {e}")
        sys.exit(1)
