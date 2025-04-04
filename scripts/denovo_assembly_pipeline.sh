# Snakemake pipeline for evaluating and refining a genome assembly without a reference genome
#ASSEMBLY = "assembly.fasta" 
NANOPORE_READS = "nanopore_reads.fq"
ILLUMINA_READS_R1 = "illumina_R1.fq"
ILLUMINA_READS_R2 = "illumina_R2.fq"

rule all:
    input:
        "assembly_stats/stats.txt",
        "busco_output/busco_summary.txt",
        "alignment/coverage_depth.txt",
        "contamination/blobtools_report.txt",
        "polished/Nextpolish_final.fasta"

# Rule 1: Assemble with shasta
rule shasta_assembler:
    input:
        reads='path_to_read/S47_trimmed.fasta'
    output:
        'path_to_output_dir/Assembly.fasta'
    shell:
        'shasta --config Nanopore-May2022 \
                --input {input.reads} \
                --command assemble \
                --assemblyDirectory {output} \
                --Reads.minReadLength 1000'

# Rule 2: Compute genome assembly statistics using QUAST
rule quast_stats:
    input:
        assembly=ASSEMBLY
    output:
        "assembly_stats/stats.txt"
    shell:
        """
        mkdir -p assembly_stats
        quast -o assembly_stats {input.assembly} > {output}
        """

# Rule 3: BUSCO completeness assessment
rule busco_assessment:
    input:
        assembly=ASSEMBLY
    output:
        "busco_output/busco_summary.txt"
    shell:
        """
        mkdir -p busco_output
        busco -i {input.assembly} -m genome -l embryophyta_odb10 -o busco_output --out_path busco_output
        """

# Rule 4: Read alignment for error estimation (Nanopore reads)
rule minimap2_alignment:
    input:
        assembly=ASSEMBLY,
        reads=NANOPORE_READS
    output:
        "alignment/aligned.sam"
    shell:
        """
        mkdir -p alignment
        minimap2 -ax map-ont {input.assembly} {input.reads} > {output}
        """

# Rule 5: Read alignment for polishing validation (Paired-End Illumina reads)
rule bwa_alignment:
    input:
        assembly=ASSEMBLY,
        reads_r1=ILLUMINA_READS_R1,
        reads_r2=ILLUMINA_READS_R2
    output:
        "alignment/illumina_aligned.bam"
    shell:
        """
        bwa index {input.assembly}
        bwa mem {input.assembly} {input.reads_r1} {input.reads_r2} | samtools sort -o {output}
        samtools index {output}
        """

# Rule 6: Compute coverage depth
rule compute_coverage:
    input:
        bam="alignment/illumina_aligned.bam"
    output:
        "alignment/coverage_depth.txt"
    shell:
        """
        samtools depth {input.bam} | awk '{{sum+=$3}} END {{print "Average depth:", sum/NR}}' > {output}
        """

# Rule 7: Contamination screening with BlobTools
rule blobtools_contamination:
    input:
        assembly=ASSEMBLY,
        reads=NANOPORE_READS
    output:
        "contamination/blobtools_report.txt"
    shell:
        """
        mkdir -p contamination
        blobtools create -i {input.assembly} -t blobtools_output -b {input.reads}
        blobtools view -i blobtools_output > {output}
        """

# Rule 8: Polishing with Racon (Long-read error correction)
rule racon_polish:
    input:
        assembly=ASSEMBLY,
        reads=NANOPORE_READS,
        alignment="alignment/aligned.sam"
    output:
        "polished/racon_polished.fasta"
    shell:
        """
        mkdir -p polished
        racon {input.reads} {input.alignment} {input.assembly} > {output}
        """

# Rule 9: Polishing with Pilon (Paired-End Illumina correction)
rule pilon_polish:
    input:
        racon_polished="polished/racon_polished.fasta",
        reads_bam="alignment/illumina_aligned.bam"
    output:
        "polished/pilon_final.fasta"
    shell:
        """
        java -Xmx16G -jar pilon.jar --genome {input.racon_polished} --fr