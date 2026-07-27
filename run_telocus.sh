#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB         
#SBATCH --time=24:00:00
#SBATCH --job-name=te_locus
#SBATCH --output=sm.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8

export TMPDIR=/scratch/user/uqctung/
export APPTAINER_TMPDIR=$TMPDIR/apptainer_tmp
export APPTAINER_CACHEDIR=$TMPDIR/apptainer_cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR

chmod +x profiles/bunya/status-sacct-robust.sh


snakemake -s telocus_pilot.smk \
    --profile profiles/bunya/ \