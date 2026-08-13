#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --job-name=mrna_te
#SBATCH --output=/scratch/user/uqctung/lmac_rna/logs/te_%j.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

# ── Locations ────────────────────────────────────────────────────────────────
QRIS_PROJECT=/QRISdata/Q9141/lmac_rna
WORKDIR=/scratch/user/uqctung/lmac_rna

mkdir -p "$WORKDIR/logs"
cd "$WORKDIR"

# ── Reference genome/annotation: symlinked from QRIS ─────────────────────────
[ -e data ]         || mkdir -p data
[ -e data/genomes ] || ln -s "$QRIS_PROJECT/data/genomes" data/genomes

# ── Bring over the already-finished fastp results (nothing on scratch yet) ──
# One-off copy-in from QRIS - this is the allowed exception under the
# Conditions of Access. Only run this once the fastp/kallisto job on QRIS
# has fully finished (check squeue / the sm.log there first).
mkdir -p results/fastp
rsync -av "$QRIS_PROJECT/results/fastp/" results/fastp/

# Sanity check nothing copied in is truncated before trusting it for STAR
for f in results/fastp/*.fastq.gz; do
    gzip -t "$f" 2>/dev/null || echo "CORRUPT FILE: $f"
done

# ── Pipeline code: copy in if not already checked out here ──────────────────
[ -e mrna_te.smk ] || cp "$QRIS_PROJECT/mrna_te.smk" .
[ -e profiles ]     || cp -r "$QRIS_PROJECT/profiles" .

export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# Containers stay on QRIS - read-only, once-off input, never written to
SIF_DIR=/QRISdata/Q9141/lmac_rna/sifs
mkdir -p "$SIF_DIR"

[ ! -f $SIF_DIR/tetools.sif ]       && apptainer pull --tmpdir "$APPTAINER_TMPDIR" $SIF_DIR/tetools.sif       docker://dfam/tetools:1.88.5
[ ! -f $SIF_DIR/star.sif ]          && apptainer pull --tmpdir "$APPTAINER_TMPDIR" $SIF_DIR/star.sif          docker://quay.io/biocontainers/star:2.7.11a--h0033a41_0
[ ! -f $SIF_DIR/samtools.sif ]      && apptainer pull --tmpdir "$APPTAINER_TMPDIR" $SIF_DIR/samtools.sif      docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0
[ ! -f $SIF_DIR/tetranscripts.sif ] && apptainer pull --tmpdir "$APPTAINER_TMPDIR" $SIF_DIR/tetranscripts.sif docker://quay.io/biocontainers/tetranscripts:2.1.3--py_0
[ ! -f $SIF_DIR/bedtools.sif ]      && apptainer pull --tmpdir "$APPTAINER_TMPDIR" $SIF_DIR/bedtools.sif      docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0

chmod +x profiles/bunya/status-sacct-robust.sh

snakemake -s mrna_te.smk --unlock --profile profiles/bunya/

snakemake -s mrna_te.smk \
    --profile profiles/bunya/ \
    --singularity-args "--bind $QRIS_PROJECT --bind $WORKDIR --bind $TMPDIR"

# ── Once-off copy of final results back to QRIS for long-term storage ───────
rsync -av --exclude 'tmp/' results/ "$QRIS_PROJECT/results/"
rsync -av lepto_batch*_gene_TE_analysis.txt "$QRIS_PROJECT/"