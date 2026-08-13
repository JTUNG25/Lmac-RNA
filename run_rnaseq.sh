#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=fastp_kallisto
#SBATCH --output=/scratch/user/uqctung/lmac_rna/logs/sm_%j.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

# ── Locations ────────────────────────────────────────────────────────────────
QRIS_PROJECT=/QRISdata/Q9141/lmac_rna
WORKDIR=/scratch/user/uqctung/lmac_rna

mkdir -p "$WORKDIR/logs"
cd "$WORKDIR"

# ── Raw inputs + reference: symlinked from QRIS, not copied ─────────────────
[ -e merged_raw ]    || ln -s "$QRIS_PROJECT/merged_raw" merged_raw
[ -e data ]          || mkdir -p data
[ -e data/genomes ]  || ln -s "$QRIS_PROJECT/data/genomes" data/genomes

# ── Pipeline code: copy in if not already checked out here ──────────────────
[ -e fastp_kal.smk ] || cp "$QRIS_PROJECT/fastp_kal.smk" .
[ -e profiles ]       || cp -r "$QRIS_PROJECT/profiles" .

# ── Reuse already-finished fastp results from QRIS ───────────────────────────
# One-off copy-in (allowed under the Conditions of Access). -u skips anything
# that already exists locally and is newer, so reruns of this script are cheap.
if [ -d "$QRIS_PROJECT/results/fastp" ]; then
    mkdir -p results/fastp
    rsync -au "$QRIS_PROJECT/results/fastp/" results/fastp/
fi

# Sanity check: make sure nothing copied in is truncated/corrupt before
# trusting it — a bad file here would otherwise silently feed kallisto garbage.
for f in results/fastp/*.fastq.gz; do
    [ -e "$f" ] || continue
    gzip -t "$f" 2>/dev/null || { echo "CORRUPT FILE, removing so it reruns: $f"; rm -f "$f"; }
done

export TMPDIR="$WORKDIR/tmp"
mkdir -p "$TMPDIR"

chmod +x profiles/bunya/status-sacct-robust.sh

# ── Run ───────────────────────────────────────────────────────────────────
snakemake -s fastp_kal.smk --profile profiles/bunya/

# ── Once-off copy of final results back to QRIS for long-term storage ───────
rsync -av --exclude 'tmp/' results/ "$QRIS_PROJECT/results/"