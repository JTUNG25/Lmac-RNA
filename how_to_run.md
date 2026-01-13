# How to Run This Workflow with Docker

This document provides instructions on how to run the RNA-seq analysis pipeline using Docker.

## 1. Installation

You will need to have Docker installed and running on your system.

## 2. Prepare Input Files

1.  **Raw FASTQ files:** Place your raw FASTQ files (e.g., `sample1.fastq.gz`) in the `data/` directory.
2.  **Reference genome:** Place your reference genome file (in FASTA format) in the `genome/` directory.
3.  **Genome annotation:** Place your genome annotation file (in GTF/GFF format) in the `genome/` directory.

## 3. Configuration

1.  **`config.yaml`:** Edit this file to specify your sample names, the paths to your reference genome and annotation files, and the contrasts for differential expression analysis.
2.  **`samples.tsv`:** Edit this file to describe your experimental design. The first column should contain the sample names (matching the `config.yaml` file), and the second column should contain the condition for each sample.

## 4. Running the Pipeline

To run the pipeline with Docker, use the following command:

```bash
snakemake --use-envmodules --cores <number_of_cores>
```

Replace `<number_of_cores>` with the number of CPU cores you want to use. Snakemake will use the `container` directives in the `Snakefile` to pull the correct Docker image for each step.