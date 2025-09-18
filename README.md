# BGCFlow (PanKB fork)

This is a fork of BGCFlow intended for use with PanKB ([repo](https://github.com/biosustain/pankb) | [website](https://pankb.org/)). For more information see [the main BGCFlow repo](https://github.com/NBChub/bgcflow) and the [BGCFlow paper](https://doi.org/10.1093/nar/gkae314).


`BGCFlow` is a systematic workflow for building pangenomes and related analytics. This fork contains additions for PanKB and a minimal workflow that runs the minimum required steps for PanKB.

At present, `BGCFlow` is only tested to work on **Linux** systems with `conda`.

## Running BGCFlow in an Azure VM

### Setting up a VM

1. **Create the VM**. Create a standard Azure x64 Linux VM with the most recent Ubuntu LTS version (general instructions [here](https://learn.microsoft.com/en-us/azure/virtual-machines/linux/quick-create-portal?tabs=ubuntu), ignore the web server section and adjust settings accordingly). During setup you can keep the VM small to save costs (you can resize CPU and RAM later). You can keep the main disk size at the default 30GB, but [`azcopy`](https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10) sometimes creates huge log files and it can be convenient to have a larger disk (though this can also be easily solved by symlinking `~/.azcopy` to somewhere in the data disk we will mount in a later step). For security and availbaility settings follow the [guidelines from DTU Biosustain](https://github.com/biosustain/guidelines) (you can also check other VMs for reference). Network settings should allow only SSH from trusted IPs (DTU IPs and any place you work remotely from).
2. **Mount a data disk**. Create and mount an additional disk for data storage. BGCFlow generates large amounts of data, so 2TB is usually the minimum to be comfortable. If you are going to be handling results from multiple runs, go for at least 4TB. Format the entire disk as `ext4`. We usually mount the data disk at `/data` (but it's not a hard requirement). Note that formatting and mounting new disks requires some general knowledge on [file systems](https://wiki.archlinux.org/title/File_systems) and the common tools for listing, partitioning, formatting and mounting disks in Linux (such as `lsblk`, `fdisk`, `mkfs.*` and the `fstab` config file).
3. **Mount the resources disk**. Mount the `bgcflow-resources` disk (from the `rg-recon` resource group). This contains mainly databases and other static files used by some of the tools of the pipeline (more details [here](https://github.com/NBChub/bgcflow/wiki/00-Installation-Guide#disk-space)). Mounting the resources disk only requires getting its UUID with `lsblk` and adding it to `fstab`. We usually mount it at `/resources` (again, not a hard requirement).
4. **Symlink the resources disk**. Clone this repository if you haven't yet (to a location of your choice). Symlink the resources directory to the mount location of the resources disk. From the root of this repository:
```bash
rm -rv ./resources
ln -sv /resources resources
```
Running `tree -L 1 ./resources` should yield something like:
```
./resources
├── automlst-simplified-wrapper-main
├── automlst-simplified-wrapper-main.back
├── checkm
├── eggnog_db
├── gtdb_download
├── gtdbtk
└── lost+found

8 directories, 0 files
```

5. **[Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)** (we usually use Miniforge). If the main VM disk is small, it is recommended you put the conda installation in the data disk rather than the default location in the home directory.
6. **Create a conda environment** with snakemake 7.31.1 and python 3.11:
```bash
conda create -n bgcflow python=3.11 bioconda::snakemake=7.31.1
```
7. **Ensure strict channel priorities for conda is disabled**:
```bash
conda config --set channel_priority disabled
```

Some additional BGCFlow installation info can be found in [the main BGCFlow wiki](https://github.com/NBChub/bgcflow/wiki/00-Installation-Guide) (note that we do **not** use the BGCFlow wrapper nor we need the GCC). 

### Running the pipeline

1. **Resize the VM**. There are no hard requirements, but a larger CPU will run faster (large runs can take a few weeks to complete). We usually use a minimum of 32 cores and 128GB of RAM for large runs.
2. **Create the config file**. Create a `config` directory and, inside it, a `config.yaml` file that should look something like this:
```yaml
taxons:
  - name: <family_name>
    reference_only: False
    # custom_samples: "config/custom_samples.csv" 

use_ncbi_data_for_checkm: True
pankb_alleleome_only: False

rules:
  classification: True
  pankb_nova: True
  pankb_minimal: True
  pankb: False
  alleleome: False
  phylons: True
```
This is also the current setting of data processing for PanKB. Replace `<family_name>` with whatever family you want to compute PanKB data for and adjust the other settings as needed. If there are custom samples to be processd, use the custom_samples property inside taxons, add the genome id as index and data stored path as path to the csv file.

3. **Create the data directory**. Create an empty directory in the data disk (any location is ok) and symlink it to `./data` (from the root of this repo), such that we effectively have an empty directory named `data` at the root of the repo. At this point, running `tree -l -L 2 .`should yield this structure:
```
.
├── CITATION.cff
├── Dockerfile
├── LICENSE
├── README.md
├── config
│   └── config.yaml
├── data -> /path/to/your/empty/data/directory
├── envs.yaml
├── resources -> /resources
│   ├── automlst-simplified-wrapper-main
│   ├── automlst-simplified-wrapper-main.back
│   ├── checkm
│   ├── eggnog_db
│   ├── gtdb_download
│   └── gtdbtk
└── workflow
    └── ...
```
You can also create the `config.yaml` file close to the data directory to keep things organized and symlink it rather than create it in the previous step.

4. **Activate the conda environment** with:
```bash
conda activate bgcflow
```

5. **Run BGCFlow** with:
```bash
snakemake --snakefile workflow/PanKB_Family --use-conda --rerun-triggers mtime -c <n_cores> --rerun-incomplete --keep-going --resources ncbi_api=1
```
Replace `<n_cores>` with the number of cores of your VM (or how many you want to allocate).

## Long-term storage of BGCFlow results

We use the `pankbpipeline` storage account in the `rg-recon` resource group to store the results and all the intermediate files from our BGCFlow runs for PanKB. Simply add it your new data as an additional directory in the runs blob and **remember to also upload the `config.yaml`**. 

You can do it easily with `azcopy` (just be wary of the huge log files issue described above). For instance, if the aggregated data directory for a specific species ran previously is `/data/<species>`, you can simply:

```bash
azcopy copy --dry-run --recursive --overwrite=false --follow-symlinks /data/<species> https://pankbpipeline.blob.core.windows.net/runs/<species>?$SAS
```

(assumes you have a `$SAS` variable with a valid [SAS token](https://learn.microsoft.com/en-us/azure/ai-services/translator/document-translation/how-to-guides/create-sas-tokens) for the Azure storage account)

With `azcopy`, alaways first do a `--dry-run` to check that everything will end up in the intended location and use `--overwrite=false` unless you specifically need to. 

## Storage of PanKB-relevant results

The files needed specifically for PanKB need to be stored the the `pankb` storage account. Instructions for that can be found in the [pankb_db README](https://github.com/biosustain/pankb_db?tab=readme-ov-file#12-copy-bgcflow-results-to-azure).
