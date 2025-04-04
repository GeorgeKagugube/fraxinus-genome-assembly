#!/usr/bin/bash -l

## State the number of cores 
#$ -pe mpi 12

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Specify the working directory here
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/scripts/assemblies/hashlr

## Activate the conda environment 
## activate the required conda environments here 
conda activate /home/uczrgwk/miniconda3/envs/haslr

## Short illumina reads
read1="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_1.fastq"
read2="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_2.fastq"

## Actual ONT data to assemble
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/S47_trimmed_fastq.fq"

# cREATE AN OUTPUT DIRECTORY HERE
output="assembly"

## Run the assembly here
haslr.py -o $output -g 800m -l $longReads -x nanopore -s $read1 $read2 -t ${NSLOTS}
