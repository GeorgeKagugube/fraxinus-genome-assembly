#!/usr/bin/bash -l

## State the cluster requirements here
#$ -pe mpi 5

# Request 12 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=12G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

## Select the run time
#$ -l h_rt=48:00:00

## Set the working directory
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/re_assemblies/uncycler_assembly/unicycler_assembly/hybrid_mode

## Activate the conda environment here 
conda activate assembly 

## Long read data to assemble
longRead="/lustre/scratch/scratch/uczrgwk/manchurica_ash/data/trimmed_DNA_ONT/S47_fastq.fq.gz"

## short read data 
read1="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_1.fastq"
read2="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_2.fastq"

## output directory
output_dir="hybrid_assembly"

## Run a hybrid assembly
unicycler -1 $read1 -2 $read2 -l $longRead -o output_dir
