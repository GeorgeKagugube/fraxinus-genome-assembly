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
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/alignments/starMapping/shasta_longReads

## Activate the conda environment
conda activate minimap2

# State the inputs here 
## Datasets 
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/combined_reads.fq.gz"
assembly="/home/uczrgwk/Scratch/manchurica_ash/assemblies/shasta_assembly/Assembly.fasta"

## Run minimap2 here 
minimap2 -ax map-ont $assembly $longReads > shasta_alignment_2.sam 

## Use samtools to sort, check coverage and basic stats from the alignment here 
samtools sort -T shasta_alignment_2.sam
samtools coverage shasta_alignment_2.sam
samtools stats shasta_alignment_2.sam
