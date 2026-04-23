#!/usr/bin/env python3

repeatmodeler = "/QRISdata/Q9141/lmac_rna/sifs/tetools.sif"
star          = "/QRISdata/Q9141/lmac_rna/sifs/star.sif"
samtools      = "/QRISdata/Q9141/lmac_rna/sifs/samtools.sif"
tetranscripts = "/QRISdata/Q9141/lmac_rna/sifs/tetranscripts.sif"
bedtools      = "/QRISdata/Q9141/lmac_rna/sifs/bedtools.sif"

import os
from pathlib import Path

(_all_r1,) = glob_wildcards("data/fastp/{sample}_R1.fastq.gz")
MRNA_SAMPLES = [s for s in _all_r1
                if os.path.exists(f"data/fastp/{s}_R2.fastq.gz")]

# Define WT vs mutant splits for TEtranscripts

WT_MRNA  = [s for s in MRNA_SAMPLES if s.startswith("WT") or s.startswith("D5")]
MUT_MRNA = [s for s in MRNA_SAMPLES if not s.startswith("WT") and not s.startswith("D5")]

GENOME    = "data/genome/JN3.fasta"
GENOME_NAME = "JN3"
GENOME_SM = "data/genome/JN3.masked.fasta"
TE_GFF    = "data/genome/JN3.te.gff3"
TE_GTF    = "data/genome/JN3.te.gtf"
TE_BED    = "data/genome/JN3.te.bed"
GENE_GTF  = "data/genome/JN3.gtf"



rule all:
    input:
        # TE annotation
        "results/repeatmodeler/JN3-families.fa",
        GENOME_SM,
        TE_GFF,
        TE_GTF,
        TE_BED,
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam",
               sample=MRNA_SAMPLES),
        "results/tetranscripts/mrna/lepto_mrna.DESeq2.TE_results.txt",
       

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — TE ANNOTATION
# ══════════════════════════════════════════════════════════════════════════════

rule repeatmodeler_build_db:
    input:
        GENOME,
    output:
        "results/repeatmodeler/JN3_db.nhr",
    params:
        prefix = "results/repeatmodeler/JN3_db",
    threads: 1
    resources:
        mem_mb  = 32000,
        runtime = 60,
    container:
        repeatmodeler
    shell:
        "BuildDatabase -name {params.prefix} -engine ncbi {input}"


rule repeatmodeler_run:
    input:
        "results/repeatmodeler/JN3_db.nhr",
    output:
        "results/repeatmodeler/JN3-families.fa",
    params:
        prefix = "results/repeatmodeler/JN3_db",
        outdir = "results/repeatmodeler",
    threads: 24
    resources:
        mem_mb  = 128000,
        runtime = 1440,
    container:
        repeatmodeler
    shell:
        """
        RepeatModeler -database {params.prefix} \
                      -threads {threads} \
                      -LTRStruct \
                      -dir {params.outdir}
        mv {params.outdir}/JN3_db-families.fa {output}
        """


rule repeatmasker:
    input:
        genome = GENOME,
        lib    = "results/repeatmodeler/JN3-families.fa",
    output:
        masked = GENOME_SM,
        gff    = "results/repeatmasker/JN3.fasta.out.gff",
    params:
        outdir  = "results/repeatmasker",
    threads: 24
    resources:
        mem_mb  = 32000,
        runtime = 480,
    container:
        repeatmodeler
    shell:
        """
        RepeatMasker -lib {input.lib} \
                     -gff -xsmall \
                     -pa {threads} \
                     -dir {params.outdir} \
                     {input.genome}
        mv {params.outdir}/JN3.fasta.masked {output.masked}
        mv {params.outdir}/JN3.fasta.out.gff {output.gff}
        """


rule make_te_gtf:
    """
    RepeatMasker GFF → GTF for TEtranscripts.
    gene_id  = TE superfamily (e.g. Gypsy, Copia, DNA/hAT …)
    transcript_id = unique locus identifier  Chr:start-end
    """
    input:
        gff = "results/repeatmasker/JN3.fasta.out.gff",
    output:
        gff = TE_GFF,
        gtf = TE_GTF,
        bed = TE_BED,
    shell:
        r"""
        # ── clean GFF3 with unique IDs ─────────────────────────────────────
        grep -v "^#" {input.gff} | awk '
        BEGIN {{ OFS="\t" }}
        {{
            # RepeatMasker GFF Target field: "Target=Motif:FamilyName start end"
            match($9, /Target=Motif:([^ ]+)/, arr)
            family = arr[1]
            if (family == "") family = "Unknown"
            locus  = $1 "_" $4 "_" $5
            print $1,$2,$3,$4,$5,$6,$7,$8,
                  "ID=" locus ";Name=" family ";family=" family
        }}' > {output.gff}

        # ── convert to TEtranscripts-compatible GTF ────────────────────────
        awk '
        BEGIN {{ OFS="\t" }}
        !/^#/ {{
            split($9, a, ";")
            locus = ""; family = ""
            for (i in a) {{
                if (a[i] ~ /^ID=/)     {{ locus  = a[i]; sub("ID=",     "", locus)  }}
                if (a[i] ~ /^family=/) {{ family = a[i]; sub("family=", "", family) }}
            }}
            print $1,$2,"exon",$4,$5,$6,$7,$8,
                  "gene_id \"" family "\"; " \
                  "transcript_id \"" locus "\"; " \
                  "family_id \"" family "\";"
        }}' {output.gff} > {output.gtf}

        # ── convert to BED for sRNA ────────────────────────
        grep -v "^#" {output.gff} | awk '
        BEGIN {{ OFS="\t" }}
        {{
            match($9, /ID=([^;]+)/,     id_arr)
            match($9, /family=([^;]+)/, fam_arr)
            locus  = id_arr[1]
            family = fam_arr[1]
            print $1, ($4 - 1), $5, locus ";" family, 0, $7
        }}' | sort -k1,1 -k2,2n > {output.bed}
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
        genome = GENOME_SM,
    output:
        index = directory("data/star_index"),
    threads: 16
    resources:
        mem_mb  = 32000,
        runtime = 60,
    container:
        star
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
      --outSAMmultNmax -1        → report ALL multi-mapping alignments
      --outFilterMultimapNmax 100 → keep reads mapping up to 100 loci
      --winAnchorMultimapNmax 200 → needed when multNmax is high
    TEtranscripts EM uses these to redistribute counts probabilistically.
    """
    input:
        r1    = "data/fastp/{sample}_R1.fastq.gz",
        r2    = "data/fastp/{sample}_R2.fastq.gz",
        index = "data/star_index",
    output:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
        log = "results/star/{sample}/Log.final.out",
    params:
        prefix = "results/star/{sample}/",
    threads: 16
    resources:
        mem_mb  = 32000,
        runtime = 120,
    container:
        star
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
             --outFileNamePrefix {params.prefix} \
             --runThreadN {threads}

        samtools index {output.bam}
        """


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — TEtranscripts DE (genes + TEs jointly)
# ══════════════════════════════════════════════════════════════════════════════

rule tetranscripts:
    """
    -t = treatment (mutant), -c = control (WT)
    --mode multi : EM redistribution of multi-mappers (required for TEs)
    --stranded   : set to 'yes'/'reverse' if stranded library, else 'no'
    """
    input:
        wt_bams  = expand("results/star/{s}/Aligned.sortedByCoord.out.bam",
                          s=WT_MRNA),
        mut_bams = expand("results/star/{s}/Aligned.sortedByCoord.out.bam",
                          s=MUT_MRNA),
        gene_gtf = GENE_GTF,
        te_gtf   = TE_GTF,
    output:
        te_res   = "results/tetranscripts/mrna/lepto_mrna.DESeq2.TE_results.txt",
        gene_res = "results/tetranscripts/mrna/lepto_mrna.DESeq2.gene_results.txt",
    params:
        outdir  = "results/tetranscripts/mrna",
        project = "lepto_mrna",
    threads: 8
    resources:
        mem_mb  = 32000,
        runtime = 240,
    container:
        tetranscripts
    shell:
        """
        TEtranscripts \
            -t {input.mut_bams} \
            -c {input.wt_bams} \
            --GTF    {input.gene_gtf} \
            --TE     {input.te_gtf} \
            --outdir {params.outdir} \
            --project {params.project} \
            --mode multi \
            --padj 0.05 \
            --minread 1 \
            --norm DESeq2 \
            --stranded no
        """