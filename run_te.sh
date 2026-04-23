#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --job-name=mrna_te
#SBATCH --output=sm.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

cd /QRISdata/Q9141/lmac_rna

# ── use /scratch (local, non-NFS) for apptainer build temp ───────────────
export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

# ── pre-pull all SIFs to NFS (persists across jobs) ──────────────────────
SIF_DIR=/QRISdata/Q9141/lmac_rna/sifs
mkdir -p $SIF_DIR

apptainer pull --tmpdir $APPTAINER_TMPDIR $SIF_DIR/tetools.sif      docker://dfam/tetools:1.88.5
apptainer pull --tmpdir $APPTAINER_TMPDIR $SIF_DIR/star.sif         docker://quay.io/biocontainers/star:2.7.11a--h0033a41_0
apptainer pull --tmpdir $APPTAINER_TMPDIR $SIF_DIR/samtools.sif     docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0
apptainer pull --tmpdir $APPTAINER_TMPDIR $SIF_DIR/tetranscripts.sif docker://quay.io/biocontainers/tetranscripts:2.1.3--py_0
apptainer pull --tmpdir $APPTAINER_TMPDIR $SIF_DIR/bedtools.sif     docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0

chmod +x profiles/bunya/status-sacct-robust.sh

snakemake -s mrna_te.smk \
    --profile profiles/bunya/ \
    --singularity-args "--bind /QRISdata/Q9141 --bind $TMPDIR"