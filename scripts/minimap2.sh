#!/usr/bin/bash -l

## State the number of cores 
#$ -pe mpi 12

## State the amount of RAM per core
#$ -l mem=10G

## Select the run time
#$ -l h_rt=48:00:00

## Specify the working directory here
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/alignments 

## Datasets 
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/trimmed_DNA_ONT/S47_trimmed.fasta"
assembly="/home/uczrgwk/Scratch/manchurica_ash/assemblies/shasta_assembly/Assembly.fasta"

## Short illumina reads
read1="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_1.fastq.gz"
read2="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_2.fastq.gz"

## Activate the conda environment
conda activate minimap2

## Run one round of Racon to polish the assembly here
minimap2 -x map-ont $assembly $longReads -t 40 > reads2asm.paf
racon $longReads reads2asm.paf $assembly -t 40 > racon_cons0.fasta

## Align raw reads to Racon polished assembly
minimap2 -ax map-ont ./racon_cons0.fasta $longReads -t 40 > reads2racon.sam
samtools view -bS -@ 40 ./reads2racon.sam -o reads2racon.bam
samtools sort -@ 40 ./reads2racon.bam -o reads2racon.sorted.bam
samtools index -@ 40 ./reads2racon.sorted.bam

