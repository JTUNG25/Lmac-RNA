#!/usr/bin/env python3

repeatmodeler = "/QRISdata/Q9141/lmac_rna/sifs/tetools.sif"
star = "/QRISdata/Q9141/lmac_rna/sifs/star.sif"
samtools = "/QRISdata/Q9141/lmac_rna/sifs/samtools.sif"
tetranscripts = "/QRISdata/Q9141/lmac_rna/sifs/tetranscripts.sif"
bedtools = "/QRISdata/Q9141/lmac_rna/sifs/bedtools.sif"

import os
from pathlib import Path

(_all_r1,) = glob_wildcards("data/fastp/{sample}_R1.fastq.gz")
MRNA_SAMPLES = [s for s in _all_r1 if os.path.exists(f"data/fastp/{s}_R2.fastq.gz")]

# ═════════════════════════
# BATCH 1: D5-5 control
# ═════════════════════════
BATCH1_CONTROL = ["D5-5-1", "D5-5-2", "D5-5-3"]
BATCH1_TREATMENT = [
    "D2-2-1",
    "D2-2-2",
    "D2-2-3",
    "D2-3-1",
    "D2-3-2",
    "D2-3-3",
    "D2-4-1",
    "D2-4-2",
    "D2-4-3",
]

# ══════════════════════════
# BATCH 2: D5-6 control
# ══════════════════════════
BATCH2_CONTROL = ["D5-6-1", "D5-6-2", "D5-6-3"]

BATCH2A_TREATMENT = [
    "A1-1-1",
    "A1-1-2",
    "A1-1-3",
    "A1-2-1",
    "A1-2-2",
    "A1-2-3",
    "A1-3-1",
    "A1-3-2",
]

BATCH2B_TREATMENT = ["A3-1", "A3-2", "A3-3"]

BATCH2C_TREATMENT = ["D1-1", "D1-2", "D1-3"]

# ══════════════════════════
# BATCH 3: D5-7 control
# ══════════════════════════
BATCH3_CONTROL = ["D5-7-1", "D5-7-2", "D5-7-3"]

BATCH3A_TREATMENT = ["R1-1", "R1-2", "R1-3"]

BATCH3B_TREATMENT = [
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
]

# ═══════════════════════════
# BATCH 4: WT-1 control
# ═══════════════════════════
BATCH4_CONTROL = ["WT-1-1", "WT-1-2", "WT-1-3"]

BATCH4A_TREATMENT = ["A13-1-1", "A13-1-2", "A13-1-3", "A13-2-1", "A13-2-2", "A13-2-3"]

BATCH4B_TREATMENT = [
    "R12-1-1",
    "R12-1-2",
    "R12-1-3",
    "R12-2-1",
    "R12-2-2",
    "R12-2-3",
    "R12-3-1",
    "R12-3-2",
    "R12-3-3",
]

GENOME = "data/genome/JN3.fasta"
GENOME_NAME = "JN3"
GENOME_SM = "data/genome/JN3.masked.fasta"
TE_GFF = "data/genome/JN3.te.gff3"
TE_GTF = "data/genome/JN3.te.gtf"
TE_BED = "data/genome/JN3.te.bed"
GENE_GTF = "data/genome/JN3.gtf"


rule all:
    input:
        "results/repeatmodeler/JN3-families.fa",
        GENOME_SM,
        TE_GFF,
        TE_GTF,
        TE_BED,
        expand(
            "results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=MRNA_SAMPLES
        ),
        expand(
            "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
            sample=MRNA_SAMPLES,
        ),
        # TEtranscripts outputs
        "results/tetranscripts/batch1/lepto_batch1_DESeq_TE_results.txt",
        "results/tetranscripts/batch2a/lepto_batch2a_DESeq_TE_results.txt",
        "results/tetranscripts/batch2b/lepto_batch2b_DESeq_TE_results.txt",
        "results/tetranscripts/batch2c/lepto_batch2c_DESeq_TE_results.txt",
        "results/tetranscripts/batch3a/lepto_batch3a_DESeq_TE_results.txt",
        "results/tetranscripts/batch3b/lepto_batch3b_DESeq_TE_results.txt",
        "results/tetranscripts/batch4a/lepto_batch4a_DESeq_TE_results.txt",
        "results/tetranscripts/batch4b/lepto_batch4b_DESeq_TE_results.txt",


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — TE ANNOTATION
# ══════════════════════════════════════════════════════════════════════════════


rule repeatmodeler_build_db:
    input:
        GENOME,
    output:
        "results/repeatmodeler/JN3_db.nhr",
    container:
        repeatmodeler
    threads: 1
    resources:
        mem_mb=32000,
        runtime=60,
    params:
        db="results/repeatmodeler/JN3_db",
    shell:
        "BuildDatabase -name {params.db} -engine ncbi {input}"


rule repeatmodeler_run:
    input:
        "results/repeatmodeler/JN3_db.nhr",
    output:
        "results/repeatmodeler/JN3-families.fa",
    container:
        repeatmodeler
    threads: 24
    resources:
        mem_mb=128000,
        runtime=1440,
    params:
        db="results/repeatmodeler/JN3_db",
        outdir="results/repeatmodeler",
    shell:
        """
        # ── run on local scratch to avoid NFS I/O bottleneck ──────────────
        SCRATCH_DIR=/scratch/temp/$SLURM_JOB_ID/repeatmodeler
        mkdir -p $SCRATCH_DIR

        # copy database files to scratch
        cp {params.db}.* $SCRATCH_DIR/

        RepeatModeler \
            -database $SCRATCH_DIR/JN3_db \
            -threads {threads} \
            -LTRStruct \
            -dir $SCRATCH_DIR

        # copy results back to NFS
        cp $SCRATCH_DIR/JN3_db-families.fa {output}
        cp $SCRATCH_DIR/consensi.fa {params.outdir}/ 2>/dev/null || true
        cp $SCRATCH_DIR/families.stk {params.outdir}/ 2>/dev/null || true
        """


rule repeatmasker_align:
    """
Run RepeatMasker alignments only (-nopost skips post-processing).
Output is JN3.fasta.cat.gz — the alignment archive.
Runs on local scratch to avoid NFS I/O bottleneck.
"""
    input:
        genome=GENOME,
        lib="results/repeatmodeler/JN3-families.fa",
    output:
        cat="results/repeatmasker/JN3.fasta.cat.gz",
        out="results/repeatmasker/JN3.fasta.out",
        gff="results/repeatmasker/JN3.fasta.out.gff",
    container:
        repeatmodeler
    threads: 24
    resources:
        mem_mb=32000,
        runtime=480,
    params:
        outdir="results/repeatmasker",
    shell:
        """
        SCRATCH_DIR=/scratch/temp/$SLURM_JOB_ID/repeatmasker
        mkdir -p $SCRATCH_DIR

        RepeatMasker \
            -lib    {input.lib} \
            -gff \
            -pa     6 \
            -dir    $SCRATCH_DIR \
            {input.genome}

        mv $SCRATCH_DIR/JN3.fasta.cat.gz {output.cat}
        mv $SCRATCH_DIR/JN3.fasta.out    {output.out}
        mv $SCRATCH_DIR/JN3.fasta.out.gff {output.gff}
        """


rule repeatmasker_postprocess:
    """
Generate soft-masked genome from existing RepeatMasker .out file.
"""
    input:
        genome=GENOME,
        out="results/repeatmasker/JN3.fasta.out",
    output:
        masked=GENOME_SM,
    container:
        bedtools
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        """
        # RepeatMasker .out starts with 3 header lines — skip them
        # columns: score div del ins query start end left strand repeat class start end left id
        tail -n +4 {input.out} | awk '{{
            print $5 "\t" ($6-1) "\t" $7
        }}' > /tmp/repeats_$SLURM_JOB_ID.bed

        bedtools maskfasta \
            -fi   {input.genome} \
            -bed  /tmp/repeats_$SLURM_JOB_ID.bed \
            -fo   {output.masked} \
            -soft

        rm /tmp/repeats_$SLURM_JOB_ID.bed
        """


rule make_te_gtf:
    """
RepeatMasker GFF → GTF for TEtranscripts.
Parses repeat name, class, family from RepeatMasker .out and .out.gff files.
TEtranscripts requires: gene_id, transcript_id, family_id, class_id
"""
    input:
        gff="results/repeatmasker/JN3.fasta.out.gff",
        out="results/repeatmasker/JN3.fasta.out",
    output:
        gff=TE_GFF,
        gtf=TE_GTF,
        bed=TE_BED,
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        r"""
        SCRATCH_DIR=/tmp/te_gtf_$$
        mkdir -p $SCRATCH_DIR

        # ── build class/family lookup from .out file ───────────────────────
        tail -n +4 {input.out} | awk '
        BEGIN {{ OFS="\t" }}
        {{
            split($11, cf, "/")
            class  = cf[1]
            family = cf[2] ? cf[2] : cf[1]
            print $10, class, family
        }}' | sort -u > $SCRATCH_DIR/te_class.txt

        # ── clean GFF3 with unique IDs ─────────────────────────────────────
        grep -v "^#" {input.gff} | awk '
        BEGIN {{ OFS="\t" }}
        {{
            match($0, /Motif:([^"]+)"/, arr)
            name  = arr[1]
            if (name == "") name = "Unknown"
            locus = $1 "_" $4 "_" $5
            print $1,$2,$3,$4,$5,$6,$7,$8,
                  "ID=" locus ";Name=" name ";family=" name
        }}' > $SCRATCH_DIR/te.gff3

        # ── convert to TEtranscripts-compatible GTF ────────────────────────
        grep -v "^#" {input.gff} | awk -v classfile=$SCRATCH_DIR/te_class.txt '
        BEGIN {{
            OFS="\t"
            while ((getline line < classfile) > 0) {{
                split(line, f, "\t")
                class_map[f[1]]  = f[2]
                family_map[f[1]] = f[3]
            }}
        }}
        {{
            match($0, /Motif:([^"]+)"/, arr)
            name   = arr[1]
            family = family_map[name] ? family_map[name] : "Unknown"
            class  = class_map[name]  ? class_map[name]  : "Unknown"
            locus  = $1 "_" $4 "_" $5
            print $1,$2,"exon",$4,$5,$6,$7,$8,
                  "gene_id \"" name "\"; " \
                  "transcript_id \"" locus "\"; " \
                  "family_id \"" family "\"; " \
                  "class_id \"" class "\";"
        }}' > $SCRATCH_DIR/te.gtf

        # ── BED for sRNA ───────────────────────────────────────────────────
        grep -v "^#" $SCRATCH_DIR/te.gff3 | awk '
        BEGIN {{ OFS="\t" }}
        {{
            match($9, /ID=([^;]+)/,     id_arr)
            match($9, /family=([^;]+)/, fam_arr)
            print $1, ($4-1), $5, id_arr[1] ";" fam_arr[1], 0, $7
        }}' | sort -k1,1 -k2,2n > $SCRATCH_DIR/te.bed

        # ── copy to NFS only at the end ────────────────────────────────────
        cp $SCRATCH_DIR/te.gff3 {output.gff}
        cp $SCRATCH_DIR/te.gtf  {output.gtf}
        cp $SCRATCH_DIR/te.bed  {output.bed}

        rm -rf $SCRATCH_DIR
        """


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — STAR INDEX + mRNA ALIGNMENT
# ══════════════════════════════════════════════════════════════════════════════


rule star_index:
    """
Build index on soft-masked genome.
genomeSAindexNbases: use 11 for ~40 Mb fungal genomes,
increase to 13 for larger assemblies (STAR will warn if wrong).
"""
    input:
        genome=GENOME_SM,
    output:
        index=directory("data/star_index"),
    container:
        star
    threads: 16
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        """
        STAR --runMode genomeGenerate \
             --genomeDir {output.index} \
             --genomeFastaFiles {input.genome} \
             --genomeSAindexNbases 11 \
             --runThreadN {threads}
        """


rule star_align_mrna:
    """
Key flags for TE analysis:
  --outSAMmultNmax -1         → report ALL multi-mapping alignments
  --outFilterMultimapNmax 100 → keep reads mapping up to 100 loci
  --winAnchorMultimapNmax 200 → needed when multNmax is high
TEtranscripts EM uses these to redistribute counts probabilistically.
"""
    input:
        r1="data/fastp/{sample}_R1.fastq.gz",
        r2="data/fastp/{sample}_R2.fastq.gz",
        index="data/star_index",
    output:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="results/star/{sample}/Log.final.out",
    container:
        star
    threads: 16
    resources:
        mem_mb=256000,
        runtime=240,
    params:
        prefix="results/star/{sample}/",
    shell:
        """
        STAR --runMode alignReads \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMstrandField intronMotif \
             --outSAMattributes NH HI AS NM \
             --outSAMmultNmax -1 \
             --outFilterMultimapNmax 100 \
             --winAnchorMultimapNmax 200 \
             --alignSoftClipAtReferenceEnds No \
             --limitBAMsortRAM 10000000000 \
             --outTmpDir /scratch/temp/$SLURM_JOB_ID/star_{wildcards.sample} \
             --outFileNamePrefix {params.prefix} \
             --runThreadN {threads}
        """


rule samtools_index:
    input:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        bai="results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
    container:
        samtools
    threads: 4
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        "samtools index -@ {threads} {input.bam}"


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — TEtranscripts DE (genes + TEs jointly)
# Separated by batch and mutant type for proper DESeq normalization
# ══════════════════════════════════════════════════════════════════════════════


# BATCH 1: D2 genes (D2-2, D2-3, D2-4) + D5-5 control
rule tetranscripts_batch1:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH1_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH1_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch1/lepto_batch1_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch1/lepto_batch1_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=480, 
    params:
        outdir="results/tetranscripts/batch1",
        project="lepto_batch1",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """


# BATCH 2A: A1 genes (A1-1, A1-2, A1-3) + D5-6 control
rule tetranscripts_batch2a:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH2A_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH2_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch2a/lepto_batch2a_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch2a/lepto_batch2a_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=480, 
    params:
        outdir="results/tetranscripts/batch2a",
        project="lepto_batch2a",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """


# BATCH 2B: A3 gene + D5-6 control
rule tetranscripts_batch2b:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH2B_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH2_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch2b/lepto_batch2b_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch2b/lepto_batch2b_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=360, 
    params:
        outdir="results/tetranscripts/batch2b",
        project="lepto_batch2b",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """


# BATCH 2C: D1 gene (DICER) + D5-6 control
rule tetranscripts_batch2c:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH2C_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH2_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch2c/lepto_batch2c_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch2c/lepto_batch2c_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=360, 
    params:
        outdir="results/tetranscripts/batch2c",
        project="lepto_batch2c",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """


# BATCH 3A: R genes R1 + D5-7 control
rule tetranscripts_batch3a:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH3A_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH3_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch3a/lepto_batch3a_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch3a/lepto_batch3a_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=720, 
    params:
        outdir="results/tetranscripts/batch3a",
        project="lepto_batch3a",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """


# BATCH 3B: R2-2, R2-3, R2-4, R2-5) + D5-7 control
rule tetranscripts_batch3b:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH3B_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH3_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch3b/lepto_batch3b_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch3b/lepto_batch3b_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=720, 
    params:
        outdir="results/tetranscripts/batch3b",
        project="lepto_batch3b",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """

# BATCH 4A: A13 genes (A13-1, A13-2) + WT-1 control
rule tetranscripts_batch4a:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH4A_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH4_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch4a/lepto_batch4a_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch4a/lepto_batch4a_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=480, 
    params:
        outdir="results/tetranscripts/batch4a",
        project="lepto_batch4a",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """


# BATCH 4B: R12 genes (R12-1, R12-2, R12-3) + WT-1 control
rule tetranscripts_batch4b:
    input:
        treatment_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH4B_TREATMENT
        ),
        control_bams=expand(
            "results/star/{s}/Aligned.sortedByCoord.out.bam", s=BATCH4_CONTROL
        ),
        gene_gtf=GENE_GTF,
        te_gtf=TE_GTF,
    output:
        te_res="results/tetranscripts/batch4b/lepto_batch4b_DESeq_TE_results.txt",
        gene_res="results/tetranscripts/batch4b/lepto_batch4b_DESeq_gene_results.txt",
    container:
        tetranscripts
    threads: 16
    resources:
        mem_mb=64000,
        runtime=720, 
    params:
        outdir="results/tetranscripts/batch4b",
        project="lepto_batch4b",
    shell:
        """
        TEtranscripts \
            -t {input.treatment_bams} \
            -c {input.control_bams} \
            --GTF {input.gene_gtf} \
            --TE {input.te_gtf} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq_default \
            --DESeq \
            --stranded no \
            --sortByPos

        mv {params.project}_DESeq_TE_results.txt   {output.te_res}
        mv {params.project}_DESeq_gene_results.txt {output.gene_res}
        """
