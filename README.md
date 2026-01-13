# RNA-seq Analysis Pipeline with Snakemake

This repository contains a Snakemake pipeline for RNA-seq data analysis. The pipeline performs quality control, trimming, alignment, quantification, and differential expression analysis. For reproducibility, the pipeline is configured to run with Singularity containers.

## 1. Installation

The most reproducible way to run this pipeline is with Singularity. You will need to have Singularity installed on your system.

The pipeline will automatically download the required container images for each step.

## 2. Prepare Input Files

1.  **Raw FASTQ files:** Place your raw FASTQ files (e.g., `sample1.fastq.gz`) in the `data/` directory.
2.  **Reference genome:** Place your reference genome file (in FASTA format) in the `genome/` directory.
3.  **Genome annotation:** Place your genome annotation file (in GTF/GFF format) in the `genome/` directory.

## 3. Configuration

1.  **`config.yaml`:** Edit this file to specify your sample names, the paths to your reference genome and annotation files, and the contrasts for differential expression analysis.
2.  **`samples.tsv`:** Edit this file to describe your experimental design. The first column should contain the sample names (matching the `config.yaml` file), and the second column should contain the condition for each sample.

## 4. Running the Pipeline

To run the pipeline with Singularity, use the following command:

```bash
snakemake --use-singularity --cores <number_of_cores>
```

Replace `<number_of_cores>` with the number of CPU cores you want to use. Snakemake will use the `container` directives in the `Snakefile` to pull the correct Singularity image for each step.

## 5. Interpreting the Results

The results of the pipeline will be saved in the `results/` directory:

*   `results/fastqc/`: Quality control reports from FastQC.
*   `results/trimmed_fastq/`: Trimmed FASTQ files.
*   `results/star_aligned/`: Aligned BAM files.
*   `results/featurecounts/`: Read counts per gene.
*   `results/deseq2/`: Differential expression analysis results in CSV format.

Each CSV file in `results/deseq2/` contains the results for a specific contrast, including log2 fold change, p-value, and adjusted p-value for each gene.
