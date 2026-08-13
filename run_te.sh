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

export PYTHONUNBUFFERED=1

# ── Locations ────────────────────────────────────────────────────────────────
QRIS_PROJECT=/QRISdata/Q9141/lmac_rna
WORKDIR=/scratch/user/uqctung/lmac_rna

mkdir -p "$WORKDIR/logs"
cd "$WORKDIR"

echo "Running from: $(pwd)"

# ── Reference genome/annotation: symlinked from QRIS ─────────────────────────
[ -e data ]         || mkdir -p data
[ -e data/genomes ] || ln -s "$QRIS_PROJECT/data/genomes" data/genomes

# ── fastp results: should already be symlinked from the earlier fastp_kal
# run (same WORKDIR). Only set it up here as a fallback if it isn't.
if [ -L results/fastp ]; then
    echo "results/fastp already symlinked, good"
elif [ -e results/fastp ]; then
    echo "results/fastp exists but is not a symlink - leaving as-is"
else
    mkdir -p results
    ln -s "$QRIS_PROJECT/results/fastp" results/fastp
fi

N_FASTP=$(ls results/fastp/*_R1.fastq.gz 2>/dev/null | wc -l)
echo "fastp R1 files visible: $N_FASTP"

# ── Pipeline code: copy in if not already checked out here ──────────────────
[ -e mrna_te.smk ] || cp "$QRIS_PROJECT/mrna_te.smk" .
[ -e profiles ]     || cp -r "$QRIS_PROJECT/profiles" .

# ── TMPDIR must be set before it's used below ────────────────────────────────
export TMPDIR="$WORKDIR/tmp"
export APPTAINER_TMPDIR="$WORKDIR/apptainer_tmp"
export APPTAINER_CACHEDIR="$WORKDIR/apptainer_cache"
mkdir -p "$TMPDIR" "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# Containers stay on QRIS - read-only, once-off input, never written to.
# Uses build --fakeroot, not pull, since pull's proot fallback fails on
# Bunya compute nodes (account not in /etc/subuid).
SIF_DIR="$QRIS_PROJECT/sifs"
mkdir -p "$SIF_DIR"

[ ! -f "$SIF_DIR/tetools.sif" ]       && apptainer build --fakeroot --tmpdir "$APPTAINER_TMPDIR" "$SIF_DIR/tetools.sif"       docker://dfam/tetools:1.88.5
[ ! -f "$SIF_DIR/star.sif" ]          && apptainer build --fakeroot --tmpdir "$APPTAINER_TMPDIR" "$SIF_DIR/star.sif"          docker://quay.io/biocontainers/star:2.7.11a--h0033a41_0
[ ! -f "$SIF_DIR/samtools.sif" ]      && apptainer build --fakeroot --tmpdir "$APPTAINER_TMPDIR" "$SIF_DIR/samtools.sif"      docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0
[ ! -f "$SIF_DIR/tetranscripts.sif" ] && apptainer build --fakeroot --tmpdir "$APPTAINER_TMPDIR" "$SIF_DIR/tetranscripts.sif" docker://quay.io/biocontainers/tetranscripts:2.1.3--py_0
[ ! -f "$SIF_DIR/bedtools.sif" ]      && apptainer build --fakeroot --tmpdir "$APPTAINER_TMPDIR" "$SIF_DIR/bedtools.sif"      docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0

chmod +x profiles/bunya/status-sacct-robust.sh

snakemake -s mrna_te.smk --unlock --profile profiles/bunya/

# --bind covers: QRIS (containers + symlinked data/genomes + symlinked
# results/fastp), WORKDIR (scratch results/logs), and TMPDIR (now correctly
# set above, no longer resolving to root).
snakemake -s mrna_te.smk \
    --profile profiles/bunya/ \
    --singularity-args "--bind $QRIS_PROJECT --bind $WORKDIR --bind $TMPDIR"

# ── Once-off copy of final results back to QRIS for long-term storage ───────
# fastp/ excluded: it's a symlink back to QRIS itself, nothing new to copy.
rsync -av --exclude 'tmp/' --exclude 'fastp/' results/ "$QRIS_PROJECT/results/"
rsync -av lepto_batch*_gene_TE_analysis.txt "$QRIS_PROJECT/" 2>/dev/null || echo "No lepto_batch*.txt files yet - run may not have completed the TEtranscripts steps"