## Code is run UCL HPC 
## Date: 30th April 2024 

## Set up the input and output directories here
input_file1="/home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/S47_trimmed_fastq.fq"
outputdir="./assembly"

## Activate the virtual environment for shasta here
conda activate /home/uczrgwk/miniconda3/envs/denovo_assembly/envs/snakemake

## Run shasta on the processed file above
## Set this up with the path to your shasta executable followed by the input and output files
shasta --config Nanopore-May2022 \
	--input $[readfile(s)] \
	--command assemble \
	--assemblyDirectory $[path_to_directory] \
	--Reads.minReadLength 1000
