## Activate the conda environment
conda activate /home/uczrgwk/miniconda3/envs/preprocessing

## Set the path to the file to be analysed here
SAMPLE="/lustre/scratch/scratch/uczrgwk/manchurica_ash/data/isoseq_data/S47_chopped.fastq"

## extract the basename for the file to use in the subequent analysis
BASE=$(basename $SAMPLE .fastq)

## Create directories to store the preprocessed qc matrices here
mkdir qc_chopped qc_chopped_trimmed qc_chopped_trimmed_filtered

## guzip the file here 
gzip $SAMPLE

fastqc -o qc_chopped $SAMPLE -t ${NSLOTS}

## Run Nanofilt here to trim the leading and trailing reads
gunzip -c ${SAMPLE}.gz | NanoFilt -q 12 --headcrop 15 --tailcrop 15 > ${BASE}_trimmed.fastq

## Assess the quality of the chopped and trimmed 
fastqc -o qc_chopped_trimmed ${BASE}_trimmed.fastq -t ${NSLOTS}

## Filter out low quality reads from using filtlong
filtlong --min_length 1000 --keep_percent 90 --target_bases 830000000 ${BASE}_trimmed.fastq > ${BASE}_trimmed_filtered.fastq

## Rerun qc on the chopped reads here
fastqc -o qc_chopped_trimmed_filtered ${BASE}_trimmed_filtered.fastq -t ${NSLOTS}
