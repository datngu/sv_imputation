#!/bin/bash
#SBATCH --job-name=sv_imputation_phase1_phase2
#SBATCH --account=nn9114k
#SBATCH --time=144:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=sv_imputation_phase1_phase2_%j.log

# Load Nextflow
module load Nextflow/24.04.2

# Set Singularity cache and tmp directories
mkdir -p $PWD/singularity
mkdir -p $PWD/singularity/tmp

export NXF_SINGULARITY_CACHEDIR="$PWD/singularity"
export SINGULARITY_CACHEDIR="$PWD/singularity"
export SINGULARITY_TMPDIR="$PWD/singularity/tmp"

echo "=== Starting Phase 1: SNP imputation ==="

nextflow run main_phase1.nf \
    -profile saga,singularity,sv_impute \
    -w work_phase1 \
    --input_manifest data/input_manifest.csv \
    --genome_build_37 data/genomes/GRCh37.fa \
    --genome_build_38 data/genomes/GRCh38.fa \
    --chain_file data/genomes/GRCh37_to_GRCh38.chain.gz \
    --snp_ref_panel_dir /cluster/projects/nn9114k/datngu/database/HC_1000G_hg38/rsid_vcfs \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir imputation_output_phase1 \
    --trace_dir imputation_traces_phase1 \
    -resume

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 1 failed. Aborting." >&2
    exit 1
fi

echo "=== Phase 1 completed successfully. Starting Phase 2: SV imputation ==="

nextflow run main_phase2.nf \
    -profile saga,singularity,sv_impute \
    -w work_phase2 \
    --phase1_snp_imputed_dir imputation_output_phase1/snp_imputed_vcfs \
    --sv_ref_panel /cluster/projects/nn9114k/datngu/database/elife_2025_sv_panel/panel.888samples.full.vcf.gz \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir imputation_output_phase2 \
    --trace_dir imputation_traces_phase2 \
    -resume

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 2 failed." >&2
    exit 1
fi

echo "=== Phase 2 completed successfully. Pipeline finished. ==="
