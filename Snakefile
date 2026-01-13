# Snakefile for RNA-seq analysis

import yaml

# Load the configuration file
configfile: "config.yaml"

# Get sample names from the config file
SAMPLES = config["samples"]

# Define container images
FASTQC_CONTAINER = "docker://biocontainers/fastqc:v0.11.9_cv8"
TRIMGALORE_CONTAINER = "docker://biocontainers/trim-galore:v0.6.6_cv1"
STAR_CONTAINER = "docker://biocontainers/star:2.7.9a--1"
SAMTOOLS_CONTAINER = "docker://biocontainers/samtools:v1.15.1_cv1"
SUBREAD_CONTAINER = "docker://biocontainers/subread:v2.0.3_cv1"
DESEQ2_CONTAINER = "docker://biocontainers/bioconductor-deseq2:v1.34.0_cv1"


# Rule all: defines the final output files of the pipeline
rule all:
    input:
        expand("results/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/trimmed_fastq/{sample}_trimmed.fq.gz", sample=SAMPLES),
        expand("results/star_aligned/{sample}.bam", sample=SAMPLES),
        expand("results/star_aligned/{sample}.bam.bai", sample=SAMPLES),
        "results/featurecounts/counts.txt",
        expand("results/deseq2/{contrast}.csv", contrast=config["contrasts"])

# Rule deseq2: perform differential expression analysis
rule deseq2:
    input:
        counts="results/featurecounts/counts.txt",
        samples_file="samples.tsv",
        script="scripts/deseq2.R"
    output:
        "results/deseq2/{contrast}.csv"
    log:
        "logs/deseq2/{contrast}.log"
    container:
        DESEQ2_CONTAINER
    shell:
        "Rscript {input.script} {input.counts} {input.samples_file} {wildcards.contrast} {output} > {log} 2>&1"

# Rule featurecounts: count reads mapped to genes
rule featurecounts:
    input:
        bams=expand("results/star_aligned/{sample}.bam", sample=SAMPLES),
        annotation=config["ref_annotation"]
    output:
        "results/featurecounts/counts.txt"
    log:
        "logs/featurecounts.log"
    container:
        SUBREAD_CONTAINER
    shell:
        "featureCounts -T 4 -a {input.annotation} -o {output} {input.bams} > {log} 2>&1"

# Rule samtools_index: index the sorted bam file
rule samtools_index:
    input:
        "results/star_aligned/{sample}.bam"
    output:
        "results/star_aligned/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    container:
        SAMTOOLS_CONTAINER
    shell:
        "samtools index {input} > {log} 2>&1"

# Rule star_align: align trimmed reads to the reference genome
rule star_align:
    input:
        fastq="results/trimmed_fastq/{sample}_trimmed.fq.gz",
        index="genome/star_index"
    output:
        "results/star_aligned/{sample}.bam"
    log:
        "logs/star_align/{sample}.log"
    container:
        STAR_CONTAINER
    shell:
        "STAR --runThreadN 4 --genomeDir {input.index} --readFilesIn {input.fastq} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results/star_aligned/{wildcards.sample} > {log} 2>&1 && mv results/star_aligned/{wildcards.sample}Aligned.sortedByCoord.out.bam {output}"

# Rule star_index: create a STAR index for the reference genome
rule star_index:
    input:
        genome=config["ref_genome"],
        annotation=config["ref_annotation"]
    output:
        directory("genome/star_index")
    log:
        "logs/star_index.log"
    container:
        STAR_CONTAINER
    shell:
        "STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.annotation} > {log} 2>&1"

# Rule trim_galore: trim adapters and low-quality reads
rule trim_galore:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/trimmed_fastq/{sample}_trimmed.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
    container:
        TRIMGALORE_CONTAINER
    shell:
        "trim_galore --fastqc -o results/trimmed_fastq {input} > {log} 2>&1"

# Rule fastqc_raw: run FastQC on raw fastq files
rule fastqc_raw:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/fastqc/{sample}_fastqc.html"
    log:
        "logs/fastqc/{sample}.log"
    container:
        FASTQC_CONTAINER
    shell:
        "fastqc {input} -o results/fastqc > {log} 2>&1"