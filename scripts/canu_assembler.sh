#!/usr/bin/bash -l

## State the number of cores 
#$ -pe mpi 12

## Request for the number of threads
#$ -pe smp 24

## State the amount of RAM per core
#$ -l mem=12G

# Set the wall time here
#$ -l h_rt=48:00:00

## Specify the working directory here
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/assemblies 

## Export path to canu2-2/bin
## Export path 
export PATH=$PATH:/home/uczrgwk/Scratch/manchurica_ash/canu-2.2/bin

## Load modules required by canu here 
module load gnuplot/6.0.rc2

## Path to long read files here
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/S47_trimmed_fastq.fq"

## Run canu assemble here
canu -d canu-assembly -p mandshurica genomeSize=803m -nanopore-raw $longReads -useGrid="false" -gridEngineResourceOption="-pe mpi THREADS -l mem=MEMORY" 
