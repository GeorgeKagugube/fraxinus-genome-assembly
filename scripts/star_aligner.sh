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
#$ -wd /home/uczrgwk/Scratch/manchurica_ash/alignments/starMapping

#You will need to install STAR using conda
module load samtools/0.1.19
module load star/2.7.3a 

## Datasets 
longReads="/home/uczrgwk/Scratch/manchurica_ash/data/trimmed_DNA_ONT/S47_trimmed.fasta"
assembly="/home/uczrgwk/Scratch/manchurica_ash/assemblies/shasta_assembly/Assembly.fasta"

## Short illumina reads
read1="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_1.fastq.gz"
read2="/home/uczrgwk/Scratch/manchurica_ash/data/short_reads/ERR4009505_2.fastq.gz"

mkdir reference

## Create an index using STAR
STAR --runMode genomeGenerate --genomeDir ./reference --genomeFastaFiles $assembly --genomeSAindexNbases 13 --runThreadN $NSLOTS;

## map the short reads to the assembly
STAR --outSAMstrandField intronMotif --genomeDir ./reference --readFilesIn $read1 $read2 --runThreadN 12 --outFileNamePrefix ./ --outSAMtype BAM SortedByCoordinate;
