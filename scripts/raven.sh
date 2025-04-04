!/usr/bin/bash -l

## State the number of cores 
#$ -pe mpi 12

## Request for the number of threads
#$ -pe smp 24

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Specify the working directory here
/home/uczrgwk/Scratch/manchurica_ash/scripts/assemblies/raven

## Set up the input and output directories here
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/S47_trimmed_fastq.fq"

## raven assembler
ravne --graphical-fragment-assembly mandchurica.gfa --threads ${NSLOTS} $longReads
