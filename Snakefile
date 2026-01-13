# Snakefile for RNA-seq analysis

import yaml
import glob
import os

# Load the configuration file
configfile: "config.yaml"

# Get sample names and paths from config
SAMPLES = config["samples"]
DATA_DIR = config["data_dir"]
SAMPLE_GROUPS = config["sample_groups"]

# Create a mapping of sample names to their fastq files
def get_fastq_files(sample):
    """Get all R1 and R2 fastq files for a sample"""
    pattern = os.path.join(DATA_DIR, f"{sample}*_R[12]*.fastq.gz")
    files = sorted(glob.glob(pattern))
    return files

def get_r1_files(sample):
    """Get all R1 fastq files for a sample"""
    pattern = os.path.join(DATA_DIR, f"{sample}*_R1*.fastq.gz")
    files = sorted(glob.glob(pattern))
    return files

def get_r2_files(sample):
    """Get all R2 fastq files for a sample"""
    pattern = os.path.join(DATA_DIR, f"{sample}*_R2*.fastq.gz")
    files = sorted(glob.glob(pattern))
    return files

# Rule all: defines the final output files of the pipeline
rule all:
    input:
        expand("results/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/trimmed_fastq/{sample}_trimmed_R1.fq.gz", sample=SAMPLES),
        expand("results/trimmed_fastq/{sample}_trimmed_R2.fq.gz", sample=SAMPLES),
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
    shell:
        "samtools index {input} > {log} 2>&1"

# Rule star_align: align trimmed reads to the reference genome
rule star_align:
    input:
        r1="results/trimmed_fastq/{sample}_trimmed_R1.fq.gz",
        r2="results/trimmed_fastq/{sample}_trimmed_R2.fq.gz",
        index="genome/star_index"
    output:
        "results/star_aligned/{sample}.bam"
    log:
        "logs/star_align/{sample}.log"
    shell:
        "STAR --runThreadN 4 --genomeDir {input.index} --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results/star_aligned/{wildcards.sample} > {log} 2>&1 && mv results/star_aligned/{wildcards.sample}Aligned.sortedByCoord.out.bam {output}"

# Rule star_index: create a STAR index for the reference genome
rule star_index:
    input:
        genome=config["ref_genome"],
        annotation=config["ref_annotation"]
    output:
        directory("genome/star_index")
    log:
        "logs/star_index.log"
    shell:
        "STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.annotation} > {log} 2>&1"

# Rule trim_galore: trim adapters and low-quality reads (paired-end)
rule trim_galore:
    input:
        r1=lambda wildcards: get_r1_files(wildcards.sample),
        r2=lambda wildcards: get_r2_files(wildcards.sample)
    output:
        r1="results/trimmed_fastq/{sample}_trimmed_R1.fq.gz",
        r2="results/trimmed_fastq/{sample}_trimmed_R2.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
    shell:
        """
        # If multiple files, concatenate them first
        if [ $(echo {input.r1} | wc -w) -gt 1 ]; then
            cat {input.r1} > /tmp/{wildcards.sample}_R1.fastq.gz
            cat {input.r2} > /tmp/{wildcards.sample}_R2.fastq.gz
            trim_galore --paired -o results/trimmed_fastq /tmp/{wildcards.sample}_R1.fastq.gz /tmp/{wildcards.sample}_R2.fastq.gz > {log} 2>&1
            rm /tmp/{wildcards.sample}_R1.fastq.gz /tmp/{wildcards.sample}_R2.fastq.gz
            # Rename output files to match expected names
            mv results/trimmed_fastq/{wildcards.sample}_R1_val_1.fq.gz {output.r1}
            mv results/trimmed_fastq/{wildcards.sample}_R2_val_2.fq.gz {output.r2}
        else
            trim_galore --paired -o results/trimmed_fastq {input.r1} {input.r2} > {log} 2>&1
            # Handle trim_galore output naming
            ls results/trimmed_fastq/{wildcards.sample}*_val_1.fq.gz | head -1 | xargs -I {{}} mv {{}} {output.r1}
            ls results/trimmed_fastq/{wildcards.sample}*_val_2.fq.gz | head -1 | xargs -I {{}} mv {{}} {output.r2}
        fi
        """

# Rule fastqc_raw: run FastQC on raw fastq files
rule fastqc_raw:
    input:
        lambda wildcards: get_fastq_files(wildcards.sample)
    output:
        "results/fastqc/{sample}_fastqc.html"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        """
        # If multiple files, concatenate them for QC
        if [ $(echo {input} | wc -w) -gt 1 ]; then
            cat {input} > /tmp/{wildcards.sample}_combined.fastq.gz
            fastqc /tmp/{wildcards.sample}_combined.fastq.gz -o results/fastqc > {log} 2>&1
            rm /tmp/{wildcards.sample}_combined.fastq.gz
            mv results/fastqc/{wildcards.sample}_combined_fastqc.html {output}
        else
            fastqc {input} -o results/fastqc > {log} 2>&1
            mv results/fastqc/*_fastqc.html {output}
        fi
        """