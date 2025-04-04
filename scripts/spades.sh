#!/usr/bin/bash -l

## State the number of cores 
#$ -pe mpi 12

## Request for the number of threads
#$ -pe smp 24

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Specify the working directory here
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/assemblies/spades_assembly

# Input files
read1="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_1.fastq.gz"
read2="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_2.fastq.gz"
long_reads="/home/uczrgwk/Scratch/manchurica_ash/data/trimmed_DNA_ONT/S47_trimmed.fasta"  # Long Oxford Nanopore reads

# Output directory for SPAdes assembly
output_dir="spades_output"

# Run SPAdes with hybrid assembly
/home/uczrgwk/Scratch/manchurica_ash/SPAdes-3.15.5-Linux/bin/spades.py --pe1-1 $read1 --pe1-2 $read2 --nanopore $long_reads -o $output_dir
