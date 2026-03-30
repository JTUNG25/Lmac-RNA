<<<<<<< HEAD
# RNA-seq Alignment Workflow for Leptosphaeria maculans (Bunya HPC)

This Snakemake workflow performs quality control and alignment of paired-end RNA-seq reads using HISAT2 on the Bunya HPC cluster using containers.

## Setup

### 1. Directory structure
Set up your project directory on Bunya:
```
/QRISdata/Q9140/lmac/lmac_rnaseq/
├── rnaseq_hisat2.smk           # Main Snakemake file
├── config.yaml                 # Configuration
├── run_rnaseq.sh              # SLURM submission script
├── profiles/bunya/             # Bunya cluster profile
│   ├── config.yaml
│   └── status-sacct-robust.sh
├── raw_data/                   # Your FASTQ files
│   ├── A1-1-1_R1.fastq.gz
│   ├── A1-1-1_R2.fastq.gz
│   └── ...
├── reference/                  # Reference genome files
│   ├── leptosphaeria_maculans_genome.fa
│   └── leptosphaeria_maculans_annotation.gtf (optional)
└── results/                    # Output directory (created automatically)
```

### 2. Configure the workflow
Edit `config.yaml` to match your file paths:
- Set `genome_fasta` to your reference genome file path
- Set `gtf_file` to your annotation file (or comment out if unavailable)

### 3. Make scripts executable
```bash
chmod +x run_rnaseq.sh
chmod +x profiles/bunya/status-sacct-robust.sh
```

## Usage

### Submit the workflow to Bunya
```bash
# Navigate to your project directory
cd /QRISdata/Q9140/lmac/lmac_rnaseq

# Submit the main job
sbatch run_rnaseq.sh
```

### Monitor progress
```bash
# Check main job status
squeue -u $USER

# Check snakemake log
tail -f snakemake_main.log

# Check individual rule logs
tail -f logs/fastqc/*.out
tail -f logs/hisat2/*.out
```

### Test the workflow (dry run)
```bash
# Load snakemake environment first
conda activate snakemake8

# Dry run to check everything
snakemake -s rnaseq_hisat2.smk --profile profiles/bunya/ -n
```

## Container Usage

This workflow uses the following containers:
- **FastQC**: `docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`
- **MultiQC**: `docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0`
- **HISAT2**: `docker://quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3`
- **Samtools**: `docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0`

All containers are automatically pulled and managed by Apptainer on Bunya.

## Resource Allocation

The workflow is optimized for Bunya with the following resource specifications:

| Rule | CPUs | Memory | Time |
|------|------|--------|------|
| FastQC | 4 | 8GB | 30min |
| Index building | 8 | 16GB | 2h |
| HISAT2 alignment | 8 | 16GB | 3h |
| BAM sorting | 4 | 8GB | 1h |
| Statistics | 2 | 4GB | 30min |
| MultiQC | 2 | 4GB | 30min |

## Output Structure

```
results/
├── fastqc/                      # FastQC quality control reports
│   ├── A1-1-1_R1_fastqc.html
│   └── ...
├── aligned/                     # Aligned BAM files
│   ├── A1-1-1.sorted.bam
│   ├── A1-1-1.sorted.bam.bai
│   └── ...
├── stats/                       # Alignment statistics
│   ├── A1-1-1.alignment_stats.txt
│   ├── alignment_summary.txt
│   └── ...
├── multiqc/                     # Aggregated QC report
│   └── multiqc_report.html
└── logs/                        # Job logs for each rule
    ├── fastqc/
    ├── hisat2/
    └── ...
```

## Key Features

1. **Automatic sample detection**: Detects all samples from FASTQ file names in `raw_data/`
2. **Container-based**: All tools run in containers for reproducibility
3. **Cluster optimized**: Designed for Bunya HPC with proper resource allocation
4. **Splice-aware alignment**: Uses GTF annotation if available
5. **Comprehensive QC**: FastQC + MultiQC reports
6. **Robust job monitoring**: Handles SLURM job status checking

## Troubleshooting

### Common issues:

1. **Permission denied for scripts**:
   ```bash
   chmod +x run_rnaseq.sh
   chmod +x profiles/bunya/status-sacct-robust.sh
   ```

2. **Container pull failures**:
   - Containers are pulled automatically; network issues may cause temporary failures
   - Jobs will retry automatically

3. **Memory issues**:
   - HISAT2 uses much less memory than STAR (16GB vs 40GB+)
   - Reduce parallel jobs if needed in `profiles/bunya/config.yaml`

4. **File not found errors**:
   - Check that `raw_data/` contains your FASTQ files
   - Verify reference genome paths in `config.yaml`

### Check job status:
```bash
# All your jobs
squeue -u $USER

# Specific job details
scontrol show job JOBID

# Check logs
tail -f logs/RULE/RULE-JOBID.out
```

## Next Steps

After alignment, you can proceed with:
- **Read counting**: Use featureCounts or HTSeq
- **Differential expression**: DESeq2, edgeR, or limma
- **Transcript assembly**: StringTie or Cufflinks
- **Quality assessment**: RSeQC or Picard tools

## Data Transfer

To download results from Bunya:
```bash
# From your local machine
scp -r USERNAME@bunya.rcc.uq.edu.au:/QRISdata/Q9140/lmac/lmac_rnaseq/results/ ./
```

## Account and Resource Notes

- Uses account: `a_qaafi_chs`
- Designed for ~74 samples (your dataset)
- Total runtime: ~8-12 hours depending on data size
- Uses Apptainer containers with GPU support enabled
- Temporary files stored in `$HOME/tmp`

=======
# Lmac-RNA
Transcriptome analasis
>>>>>>> ecd37a6f4db7277780324cd06a9182d5d535b276
