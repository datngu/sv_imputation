#!/bin/bash
#SBATCH --job-name=sv_imputation_phase2
#SBATCH --account=nn9114k
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=sv_imputation_phase2_%j.log

# Load Nextflow
module load Nextflow/24.04.2

# Set Singularity cache and tmp directories
mkdir -p /cluster/projects/nn9114k/datngu/singularity/tmp

export SINGULARITY_CACHEDIR="/cluster/projects/nn9114k/datngu/singularity"
export SINGULARITY_TMPDIR="/cluster/projects/nn9114k/datngu/singularity/tmp"
export NXF_SINGULARITY_CACHEDIR="/cluster/projects/nn9114k/datngu/singularity"

nextflow run main_phase2.nf \
    -profile saga,singularity,sv_impute \
    --phase1_snp_imputed_dir imputation_output_phase1/snp_imputed_vcfs \
    --sv_ref_panel /cluster/projects/nn9114k/datngu/database/elife_2025_sv_panel/panel.888samples.full.vcf.gz \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir imputation_output_phase2 \
    --trace_dir imputation_traces_phase2 \
    -resume
