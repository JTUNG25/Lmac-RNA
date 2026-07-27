# ══════════════════════════════════════════════════════════════════════════
# PILOT: locus-resolved TE quantification via unique-mapping reads
#
# Does NOT require a new alignment or new GTF. Reuses:
#   - existing STAR BAMs (already tagged with NH per read)
#   - existing TE_GTF (already has transcript_id = "chr_start_end" per locus,
#     alongside gene_id = subfamily name, from the make_te_gtf rule)
#
# Logic: keep only NH:i:1 (uniquely-mapping) reads, then count them against
# the TE GTF grouped by transcript_id (locus) instead of gene_id (subfamily).
# Loci too similar to their neighbours will come back near-zero — that means
# "unresolvable with unique reads," not "unexpressed."
# ══════════════════════════════════════════════════════════════════════════

featurecounts = "/QRISdata/Q9141/lmac_rna/sifs/subread.sif"

samtools = "/QRISdata/Q9141/lmac_rna/sifs/samtools.sif"

# ── Pilot scope: change these two lines to pilot a different family/batch ──
PILOT_FAMILY = "Penelope"
PILOT_SAMPLES = [
    "D5-6-1", "D5-6-2", "D5-6-3",
    "A1-1-1", "A1-1-2", "A1-1-3",
    "A1-2-1", "A1-2-2", "A1-2-3",
    "A1-3-1", "A1-3-2",
]  # Δago1 + its D5-6 control
GENOME = "data/genome/JN3.fasta"
GENOME_NAME = "JN3"
GENOME_SM = "data/genome/JN3.masked.fasta"
TE_GFF = "data/genome/JN3.te.gff3"
TE_GTF = "data/genome/JN3.te.gtf"
TE_BED = "data/genome/JN3.te.bed"
GENE_GTF = "data/genome/JN3.gtf"
# ─────────────────────────────────────────────────────────────────────────

TE_GTF_PILOT = f"data/genome/JN3.te.{PILOT_FAMILY}.gtf"


rule all:
    input:
        f"results/featurecounts_pilot/{PILOT_FAMILY}_locus_counts.txt",
        f"results/featurecounts_pilot/{PILOT_FAMILY}_subfamily_counts.txt",


rule filter_te_gtf_pilot_family:
    """
    Subset the full TE GTF down to one family, so the pilot count matrix is
    small enough to sanity-check by eye before considering a genome-wide run.
    """
    input:
        gtf=TE_GTF,
    output:
        gtf=TE_GTF_PILOT,
    resources:
        mem_mb=2000,
        runtime=10,
    shell:
        """
        awk -F'\t' -v fam='family_id "{PILOT_FAMILY}"' \
            '$9 ~ fam' {input.gtf} > {output.gtf}
        """


rule filter_unique_reads:
    """
    Keep only uniquely-mapping reads (NH:i:1) from the existing multi-mapped
    STAR BAMs. No realignment needed — STAR already tags NH per read from the
    original star_align_mrna run.

    If samtools < 1.10 (no -e filter expression support), replace the
    samtools view line with:
        samtools view -h {input.bam} \
          | awk '$0 ~ /^@/ || $0 ~ /NH:i:1\\b/' \
          | samtools view -b -o {output.bam} -
    """
    input:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        bam="results/star_unique/{sample}/Aligned.unique.bam",
        bai="results/star_unique/{sample}/Aligned.unique.bam.bai",
    container:
        samtools
    threads: 4
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view -@ {threads} -b -e '[NH]==1' {input.bam} > {output.bam}
        samtools index -@ {threads} {output.bam}
        """


rule featurecounts_locus_pilot:
    """
    Locus-resolved counts using unique-mapping reads only.
    -g transcript_id groups by individual genomic locus (chr_start_end) —
    the same GTF you already use for TEtranscripts/degeneration analysis,
    just grouped by a different attribute.
    """
    input:
        bams=expand("results/star_unique/{s}/Aligned.unique.bam", s=PILOT_SAMPLES),
        gtf=TE_GTF_PILOT,
    output:
        counts=f"results/featurecounts_pilot/{PILOT_FAMILY}_locus_counts.txt",
        summary=f"results/featurecounts_pilot/{PILOT_FAMILY}_locus_counts.txt.summary",
    container:
        featurecounts
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
    shell:
        """
        mkdir -p $(dirname {output.counts})
        featureCounts \
            -a {input.gtf} \
            -o {output.counts} \
            -g transcript_id \
            -t exon \
            -p \
            --countReadPairs \
            -T {threads} \
            {input.bams}
        """
        # Note: if your featureCounts version is >= 2.0.3, add --countReadPairs
        # alongside -p for the current fragment-counting behaviour; check
        # `featureCounts -v` first since the flag doesn't exist in older builds.


rule featurecounts_subfamily_pilot:
    """
    Same reads, same GTF, but grouped by gene_id (subfamily) instead — a
    same-input comparison point against your existing TEtranscripts subfamily
    counts. Won't match exactly (this is unique-reads-only + featureCounts'
    counting logic, not TEtranscripts' EM redistribution), but useful as a
    sanity check on how much signal survives the unique-only filter overall.
    """
    input:
        bams=expand("results/star_unique/{s}/Aligned.unique.bam", s=PILOT_SAMPLES),
        gtf=TE_GTF_PILOT,
    output:
        counts=f"results/featurecounts_pilot/{PILOT_FAMILY}_subfamily_counts.txt",
    container:
        featurecounts
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
    shell:
        """
        mkdir -p $(dirname {output.counts})
        featureCounts \
            -a {input.gtf} \
            -o {output.counts} \
            -g gene_id \
            -t exon \
            -p \
            --countReadPairs \
            -T {threads} \
            {input.bams}
        """

