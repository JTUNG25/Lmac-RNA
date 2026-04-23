#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=mrna_te
#SBATCH --output=sm.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

# Set working directory (change this to your project directory)
cd /QRISdata/Q9141/lmac_rna

export TMPDIR=/QRISdata/Q9141/lmac_rna/tmp
mkdir -p $TMPDIR

# Make status script executable
chmod +x profiles/bunya/status-sacct-robust.sh

# Run snakemake
snakemake -s rnaseq_hisat2.smk --profile profiles/bunya/ 
