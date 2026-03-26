#!/bin/bash
#SBATCH --job-name=sv_imputation_1stage
#SBATCH --account=nn9114k
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=sv_imputation_1stage_%j.log

# Load Nextflow
module load Nextflow/24.04.2

# Set Singularity cache and tmp directories
mkdir -p $PWD/singularity
mkdir -p $PWD/singularity/tmp

export NXF_SINGULARITY_CACHEDIR="$PWD/singularity"
export SINGULARITY_CACHEDIR="$PWD/singularity"
export SINGULARITY_TMPDIR="$PWD/singularity/tmp"

nextflow run main_1_stage.nf \
    -profile saga,singularity,sv_impute \
    --input_manifest data/input_manifest.csv \
    --genome_build_37 data/genomes/GRCh37.fa \
    --genome_build_38 data/genomes/GRCh38.fa \
    --chain_file data/genomes/GRCh37_to_GRCh38.chain.gz \
    --sv_ref_panel /cluster/projects/nn9114k/datngu/database/elife_2025_sv_panel/panel.888samples.full.vcf.gz \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir imputation_output_1stage \
    --trace_dir imputation_traces_1stage \
    -resume
