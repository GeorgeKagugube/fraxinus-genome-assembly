#!/usr/bin/bash -l

## State the number of cores 
#$ -pe mpi 12

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Specify the working directory here
#$ -wd /lustre/scratch/scratch/uczrgwk/manchurica_ash/scripts/assemblies/busco_analysis 

## activate the required conda environments here 
#conda activate /home/uczrgwk/miniconda3/envs/busco

assembly="/home/uczrgwk/Scratch/manchurica_ash/assemblies/flye_assembly/assembly.fasta"

## Run busco here
/home/uczrgwk/Scratch/manchurica_ash/software/busco-4.0.1/bin/busco -i $assembly -o qc_analysis -l eukaryota_odb10 -m geno -f --scaffold_composition
/home/uczrgwk/Scratch/manchurica_ash/software/busco-4.0.1/bin/busco -i $assembly -o qc_analysis -l eudicots_odb10 -m geno -f --scaffold_composition
