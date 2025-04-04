#!/usr/bin/bash -l
## State the number of cores 
#$ -pe mpi 16

## Request for the number of threads
#$ -pe smp 24

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Specify the working directory here
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/scripts/assemblies/redbeans 

## Long read data to assemble
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/isoseq_data/S47_trimmed_fastq.fq"

## Run redbeans here
/home/uczrgwk/Scratch/manchurica_ash/software/wtdbg2/wtdbg2 -x ont -g 800m -t ${NSLOTS} -i $longReads -fo mandchurica 
#/home/uczrgwk/Scratch/manchurica_ash/software/wtdbg2/wtdbg2.pl -o manchurica -g 800000000 -t ${NSLOTS} -x ont $longReads
#/home/uczrgwk/Scratch/manchurica_ash/software/wtdbg2/wtdbg2 -x ont -g 800000000 -t ${NSLOTS} -i $longReads -fo mandchurica
/home/uczrgwk/Scratch/manchurica_ash/software/wtdbg2/wtpoa-cns -t ${NSLOTS} -i mandchurica.ctg.lay.gz -fo mandchurica.ctg.fa
