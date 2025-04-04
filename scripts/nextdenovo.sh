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
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/scripts/assemblies/nextDenovo/NextDenovo

# activate teh conda environment here 
conda activate /home/uczrgwk/miniconda3/envs/nextdenovo

## Set the input and output files here
longReads="/home/uczrgwk/Scratch/manchurica_ash/scripts/assemblies/nextDenovo/NextDenovo/S47_trimmed_fastq.fq"

## Prepare the input file to nextDenovo
ls $longReads > ./input.fofn

##
cp doc/run.cfg ./

## 
nextDenovo genome_size=800000000 read_type=ont run.cfg 
