# ECHO Genome Annotation

ECHO computes comprehensive genome annotations by combining and filtering the results of EVM/PASA, GETA, Mikado, and Helixer. 
ECHO is designed to run via the SLURM job scheduler on a cloud or other high-performance computing system.

ECHO runs in two phases, the annotation phase and the filtering phase. 
The annotation phase processes all the provided sources of evidence (RNAseq, protein databases, neighbor species, and *ab initio* preditcions), computes multiple draft genome annotations, and merges them into a single draft. 
The filtering phase uses a semi-supervised random forest algorithm to remove spurious gene models.

The pipline also includes **TidyGFF**, a robust script for formatting GFF files accoding to user defined naming conventions.
After running the pipline, use ```bin/TidyGFF.py``` to ready your GFF for public distribution.

## Installation and environment
---------------
The quickest way to obtain ECHO is to download the singularity image, clone the repo, and install Snakemake with conda.
```bash
# Download the singularity image.
wget ...

# Download the repo
git clone git@github.com:plantgenomicslab/ECHO.git
cd ECHO

# Create an environment with Snakemake version 7
conda create -n snake -c bioconda -c conda-forge \
  python=3.11 \
  snakemake=7 -y

conda activate snake
```
If you have root permission and you prefer to build the singularity image from source:
```bash
# Build the singularity image
cd ECHO/singularity

sudo singularity build \
  --bind ~/ECHO/singularity/.bashrc:/.bashrc \
  echo.sif echo.def
```

## Preparing for a run
-----------------

### Input data download
1. **Download RNAseq data (fastq - paired/unpaired/both):** 
All fastq files (all library types) should be placed into a single folder named with '.fastq.gz' suffix. Paired-end reads should be named '_1.fastq.gz' and '_2.fastq.gz'. Symlinks may be used to avoid storing multiple copies of the raw data. 

2. **Download proteins (fasta):**
Download protein fasta files from Uniprot, OrthoDB, etc.

3. **Download neighbor species annotations (fasta, gff3):** 
Place all assembly (fasta) files into one folder and all annotation files (GFF3) into another. I like to organize these into ```neighbors/gff``` and ```neighbors/fasta```. If possible, obtain the CDS sequences as well for use with EDTA, otherwise you will need to generate these yourself.

4. **Download RexDB:**
In addition to EDTA output, the filter requires a library of plant-specific known repeat elements from RexDB.
```bash
wget https://github.com/repeatexplorer/rexdb/blob/main/Viridiplantae_v4.0.fasta
```

### Construct a repeat library

We use EDTA to construct a repeat sequence annotation for the genome. For EDTA, obtain cds sequences for the neighbor species you will use in the analysis. Then it can be run via the singularity image:
```
singularity exec [echo.sif] EDTA.pl \
	--genome [genome.fasta] \
	--cds [neighbors.cds] \
	--anno 1 \
	--threads [threads]
```

### Configure the annotation and filter phases
Both the annotation and filtering phases are managed by YAML control files (see ```config``` folder for examples). Customize the control files by first copy them from ```config``` to the ```ECHO``` folder:
```bash
cp config/config_annotate.yml config/config_filter.yml .
```
Then fill in the necessary fields as appropriate. Aside from providing input data, you will need to add your SLURM account information.


## Running the pipeline
-----------------
Before running each phase, it is a good idea to perform a dry run with snakemake to ensure that the configuration was done properly and all file paths are correct:
```
snakemake -np -s bin/Snakefile_annotate
```

Currently, the tool is only set up to work with the SLURM job scheduler. Run the annotation phase like so:
```bash
# Run the annotation phase
sbatch -A [account] \
    -p [partition] \
    -c 1 \
    --mem=1g \
    -J annotate \
    -o annotate.out \
    -e annotate.err \
    --wrap="bin/annotate.sh"
```
The final output of the annotation phase is called ```complete_draft.gff3```. This is represents the combined gff file from EVM/PASA, Helixer, Mikado, and GETA. Filtration can procede from this combined output or any of the intermediate draft annotations.

Check the filter configuration by performing a dry run and then execute ```bin/filter.sh``` via SLURM.

```bash
# Dry run for the filter phase
snakemake -np -s bin/Snakefile_filter

# Run the filter phase
sbatch -A [account] \
    -p [partition] \
    -c 1 \
    --mem=4g \
    -J filter \
    -o filter.out \
    -e filter.err \
    --wrap="bin/filter.sh"
```

The final output of the filter is ```ECHO/FILTER/filter.gff3```.

### Tuning the filter
After the filtration script has completed, the results may be fine-tuned by adjusting the parameters of the random forest. 
To adjust hyperparameters, rerun the python script outside of Snakemake:
```
singularity exec [echo.sif] python bin/Filter.py --help
``` 

### Logs and troubleshooting 

Given variations in input data size, the program may terminate prematurely due to memory issues. In the event that the program stops prematurely, first check the SLURM log to identify which Snakemake rule produced the error.
```
grep Error [annotation.err|filter.err]
```
Log files for each rule are stored in ```ECHO/logs```. Look in ```ECHO/logs/[rule].err``` and/or ```ECHO/logs/[rule].out``` for potential memory issues. If a rule is crashing due to memory constraints, the memory can be adjusted in the  ```Cluster Configuration``` section of each configuration file.

## Formatting output with TidyGFF

```
singularity exec [echo.sif] python bin/TidyGFF.py \
    [pre], \
    [gff], \
    --out [output_name] \
    --splice-name [t|mRNA|etc.] \
    --justify [digit_number] \
    --sort \
    --chrom-regex [chromosome_regex] \
    --source [source]
```
