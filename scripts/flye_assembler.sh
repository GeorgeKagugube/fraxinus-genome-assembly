#!/usr/bin/bash -l

## Ask for grid resources here 
## State the number of cores 
#$ -pe mpi 16

## Request for the number of threads
#$ -pe smp 24

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Specify the working directory here
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/scripts/assemblies/flye

# activate teh conda environment here 
conda activate /home/uczrgwk/miniconda3/envs/flye_assembly

## Set the input and output files here
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/S47_trimmed_fastq.fq"

# Output directory
outdir="flye"

## Run flye here
flye --nano-raw $longReads -g 800m -o $outdir -t ${NSLOTS} --meta --keep-haplotypes --scaffold 

