#!/usr/bin/env python3

# Container definitions
fastqc_container = "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
multiqc_container = "docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
hisat2_container = "docker://quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3"
samtools_container = "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"

GENOME = "data/genome/JN3.fasta"
GTF = "data/genome/JN3.gtf"


(all_samples,) = glob_wildcards("concatenated_fastq/{sample}_R1.fastq.gz")
SAMPLES = all_samples


rule target:
    input:
        # Quality control
        expand(
            "results/fastqc/{sample}_R{read}_fastqc.html", sample=SAMPLES, read=[1, 2]
        ),
        "results/multiqc/multiqc_report.html",
        # Alignment
        expand("results/hisat2/{sample}.sorted.bam", sample=SAMPLES),
        expand("results/hisat2/{sample}.sorted.bam.bai", sample=SAMPLES),
        # Alignment statistics
        expand("results/stats/{sample}.alignment_stats.txt", sample=SAMPLES),
        "results/stats/alignment_summary.txt",


rule fastqc:
    input:
        r1="concatenated_fastq/{sample}_R1.fastq.gz",
        r2="concatenated_fastq/{sample}_R2.fastq.gz",
    output:
        r1_html="results/fastqc/{sample}_R1_fastqc.html",
        r1_zip="results/fastqc/{sample}_R1_fastqc.zip",
        r2_html="results/fastqc/{sample}_R2_fastqc.html",
        r2_zip="results/fastqc/{sample}_R2_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    container:
        fastqc_container
    threads: 4
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        outdir="results/fastqc",
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2}
        """

# Optional: Build HISAT2 index with splice sites (if GTF available)
rule build_hisat2_index_with_ss:
    input:
        genome=GENOME,
        gtf=GTF,
    output:
        ss="reference/hisat2_index/splicesites.txt",
        exon="reference/hisat2_index/exons.txt",
        idx=expand("reference/hisat2_index/genome_ss.{i}.ht2", i=range(1, 9)),
    log:
        "logs/hisat2_index/build_index_ss.log",
    container:
        hisat2_container
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        prefix="reference/hisat2_index/genome_ss",
    shell:
        """
        mkdir -p reference/hisat2_index
        
        # Extract splice sites and exons
        hisat2_extract_splice_sites.py {input.gtf} > {output.ss}
        hisat2_extract_exons.py {input.gtf} > {output.exon}
        
        # Build index with splice site information
        hisat2-build -p {threads} \
            --ss {output.ss} \
            --exon {output.exon} \
            {input.genome} {params.prefix}
        """


rule hisat2_align:
    input:
        r1="concatenated_fastq/{sample}_R1.fastq.gz",
        r2="concatenated_fastq/{sample}_R2.fastq.gz",
        idx=expand("reference/hisat2_index/genome_ss.{i}.ht2", i=range(1, 9)),
    output:
        sam=temp("results/hisat2/{sample}.sam"),
        summary="results/hisat2/{sample}.hisat2_summary.txt",
    log:
        "logs/hisat2/{sample}.log",
    container:
        hisat2_container
    threads: 8
    resources:
        mem_mb=32000,
        runtime=180,
    params:
        index_prefix="reference/hisat2_index/genome_ss",
        rg_id="{sample}",
        rg_sm="{sample}",
    shell:
        """
        hisat2 -p {threads} \
            -x {params.index_prefix} \
            -1 {input.r1} -2 {input.r2} \
            --rg-id {params.rg_id} --rg SM:{params.rg_sm} \
            --summary-file {output.summary} \
            -S {output.sam}
        """


rule sam_to_bam:
    input:
        "results/hisat2/{sample}.sam",
    output:
        "results/hisat2/{sample}.bam",
    log:
        "logs/sam_to_bam/{sample}.log",
    container:
        samtools_container
    threads: 2
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        "samtools view -bS {input} > {output}"


rule sort_bam:
    input:
        "results/hisat2/{sample}.bam",
    output:
        bam="results/hisat2/{sample}.sorted.bam",
        bai="results/hisat2/{sample}.sorted.bam.bai",
    log:
        "logs/samtools_sort/{sample}.log",
    container:
        samtools_container
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """


# Generate alignment statistics
rule alignment_stats:
    input:
        "results/hisat2/{sample}.sorted.bam",
    output:
        stats="results/stats/{sample}.alignment_stats.txt",
        flagstat="results/stats/{sample}.flagstat.txt",
        idxstats="results/stats/{sample}.idxstats.txt",
    log:
        "logs/samtools_stats/{sample}.log",
    container:
        samtools_container
    threads: 2
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        """
        mkdir -p results/stats
        
        # Basic statistics
        samtools stats {input} > {output.stats}
        
        # Flag statistics
        samtools flagstat {input} > {output.flagstat}
        
        # Index statistics
        samtools idxstats {input} > {output.idxstats}
        """


# Aggregate alignment statistics
rule aggregate_stats:
    input:
        expand("results/stats/{sample}.flagstat.txt", sample=SAMPLES),
    output:
        "results/stats/alignment_summary.txt",
    log:
        "logs/aggregate_stats/aggregate.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    run:
        with open(output[0], "w") as out:
            out.write(
                "Sample\tTotal_reads\tMapped_reads\tMapping_rate\tProperly_paired\n"
            )

            for sample in SAMPLES:
                flagstat_file = f"results/stats/{sample}.flagstat.txt"

                with open(flagstat_file, "r") as f:
                    lines = f.readlines()

                # Parse flagstat output
                total_reads = lines[0].split()[0]
                mapped_line = [l for l in lines if "mapped (" in l][0]
                mapped_reads = mapped_line.split()[0]
                mapping_rate = mapped_line.split("(")[1].split("%")[0]

                properly_paired_line = [l for l in lines if "properly paired" in l][
                    0
                ]
                properly_paired = properly_paired_line.split("(")[1].split("%")[0]

                out.write(
                    f"{sample}\t{total_reads}\t{mapped_reads}\t{mapping_rate}%\t{properly_paired}%\n"
                )


# MultiQC report
rule multiqc:
    input:
        fastqc=expand(
            "results/fastqc/{sample}_R{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]
        ),
        hisat2=expand("results/hisat2/{sample}.hisat2_summary.txt", sample=SAMPLES),
        stats=expand("results/stats/{sample}.alignment_stats.txt", sample=SAMPLES),
    output:
        "results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc/multiqc.log",
    container:
        multiqc_container
    threads: 2
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        outdir="results/multiqc",
        indir="results",
    shell:
        """
        mkdir -p {params.outdir}
        multiqc -o {params.outdir} {params.indir}
        """


# Clean up intermediate files
rule clean:
    shell:
        """
        rm -f results/hisat2/*.bam
        echo "Cleaned unsorted BAM files"
        """
