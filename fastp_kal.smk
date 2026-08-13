#!/usr/bin/env python3


fastp = "docker://quay.io/biocontainers/fastp:1.0.1--heae3180_0"
kallisto = "docker://quay.io/biocontainers/kallisto:0.51.1--h2b92561_2"
multiqc = "docker://quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0"

# TODO: confirm this points at the NEW genome assembly / annotation-derived
# transcriptome fasta, not a placeholder left over from the JN3 build.
TRANSCRIPTOME = "data/genomes/Lmac_D5_transcripts.fasta"

SAMPLES = [
    # ── New samples — batch1 (control: D5-5) ────────────────────────────────
    "D2-2-1",
    "D2-2-2",
    "D2-2-3",
    "D2-3-1",
    "D2-3-2",
    "D2-3-3",
    "D5-5-1",
    "D5-5-2",
    "D5-5-3",
    # ── New samples — batch2 (control: D5-6) ────────────────────────────────
    "A1-1-1",
    "A1-1-2",
    "A1-1-3",
    "A1-2-1",
    "A1-2-2",
    "A1-2-3",
    "A3-1",
    "A3-2",
    "A3-3",
    "D1-1",
    "D1-2",
    "D1-3",
    "D5-6-1",
    "D5-6-2",
    "D5-6-3",
    # ── New samples — batch3 (control: D5-7) ────────────────────────────────
    "R1-1",
    "R1-2",
    "R1-3",
    "R2-2-1",
    "R2-2-2",
    "R2-2-3",
    "R2-3-1",
    "R2-3-2",
    "R2-3-3",
    "R2-4-1",
    "R2-4-2",
    "R2-4-3",
    "R2-5-1",
    "R2-5-2",
    "R2-5-3",
    "D5-7-1",
    "D5-7-2",
    "D5-7-3",
    # ── New samples — batch4 (control: WT-1) ────────────────────────────────
    "A13-1-1",
    "A13-1-2",
    "A13-1-3",
    "A13-2-1",
    "A13-2-2",
    "A13-2-3",
    "R12-1-1",
    "R12-1-2",
    "R12-1-3",
    "R12-2-1",
    "R12-2-2",
    "R12-2-3",
    "R12-3-1",
    "R12-3-2",
    "R12-3-3",
    "WT-1-1",
    "WT-1-2",
    "WT-1-3",
    # (batch5)
    "A2-1-1",
    "A2-1-2",
    "A2-1-3",
    "A2-2-1",
    "A2-2-2",
    "A2-2-3",
    "A2-3-1",
    "A2-3-2",
    "A2-3-3",
    "R3-1-1",
    "R3-1-2",
    "R3-1-3",
    "R3-2-1",
    "R3-2-2",
    "R3-2-3",
    "R3-3-1",
    "R3-3-2",
    "R3-3-3",
    "R1-2-1",
    "R1-2-2",
    "R1-2-3",
    "WT-2-1",
    "WT-2-2",
    "WT-2-3",
]


rule all:
    input:
        expand("results/fastp/{sample}.html", sample=SAMPLES),
        expand("results/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
        "results/multiqc/multiqc_report.html",


# ── Build Kallisto index ─────────────────────────────────────────────────────
rule kallisto_index:
    input:
        TRANSCRIPTOME,
    output:
        "reference/kallisto/transcriptome.idx",
    log:
        "logs/kallisto_index/build.log",
    container:
        kallisto
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        """
        mkdir -p reference/kallisto
        kallisto index \
            -i {output} \
            {input} \
            2>{log}
        """


# ── Adapter trimming and QC with fastp ───────────────────────────────────────
rule fastp:
    input:
        r1="merged_raw/{sample}_R1.fastq.gz",
        r2="merged_raw/{sample}_R2.fastq.gz",
    output:
        r1="results/fastp/{sample}_R1.fastq.gz",
        r2="results/fastp/{sample}_R2.fastq.gz",
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    container:
        fastp
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        """
        mkdir -p results/fastp
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --length_required 35 \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 20 \
            2>{log}
        """


rule kallisto_quant:
    input:
        r1="results/fastp/{sample}_R1.fastq.gz",
        r2="results/fastp/{sample}_R2.fastq.gz",
        idx="reference/kallisto/transcriptome.idx",
    output:
        abundance="results/kallisto/{sample}/abundance.tsv",
        h5="results/kallisto/{sample}/abundance.h5",
        run_info="results/kallisto/{sample}/run_info.json",
    log:
        "logs/kallisto/{sample}.log",
    container:
        kallisto
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        outdir="results/kallisto/{sample}",
        # bootstrap samples for uncertainty estimation — useful for sleuth,
        # harmless for tximport/DESeq2; remove -b if you want faster runs
        bootstrap=100,
    shell:
        """
        mkdir -p {params.outdir}
        kallisto quant \
            -i {input.idx} \
            -o {params.outdir} \
            -b {params.bootstrap} \
            --threads {threads} \
            {input.r1} {input.r2} \
            2>{log}
        """


rule multiqc:
    input:
        fastp_json=expand("results/fastp/{sample}.json", sample=SAMPLES),
        kallisto_log=expand("logs/kallisto/{sample}.log", sample=SAMPLES),
    output:
        "results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc/multiqc.log",
    container:
        multiqc
    threads: 1
    resources:
        mem_mb=4000,
        runtime=20,
    params:
        outdir="results/multiqc",
        search_dirs="results/fastp logs/kallisto",
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p tmp
        TMPDIR=$(pwd)/tmp \
            multiqc \
            {params.search_dirs} \
            --outdir {params.outdir} \
            --force \
            2>{log}
        """
