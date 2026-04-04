#!/bin/bash
#SBATCH --job-name=sv_imputation_1stage
#SBATCH --account=p33_norment
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=sv_imputation_1stage_%j.log

# Load Nextflow
module load Nextflow/24.04.2

# Set Singularity cache to pre-downloaded container directory (no internet on TSD)
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/p33/users/datn/containers
export SINGULARITY_CACHEDIR=/cluster/projects/p33/users/datn/containers
export SINGULARITY_TMPDIR=/cluster/projects/p33/users/datn/containers/tmp

mkdir -p $SINGULARITY_TMPDIR

nextflow run main_1_stage.nf \
    -profile tsd,singularity,sv_impute \
    -w work_1stage \
    --input_manifest data/input_manifest.csv \
    --genome_build_37 /cluster/projects/p33/users/datn/data/genomes/GRCh37.fa \
    --genome_build_38 /cluster/projects/p33/users/datn/data/genomes/GRCh38.fa \
    --chain_file /cluster/projects/p33/users/datn/data/genomes/GRCh37_to_GRCh38.chain.gz \
    --sv_ref_panel /cluster/projects/p33/users/datn/data/elife_2025_sv_panel/panel.888samples.full.vcf.gz \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir imputation_output_1stage \
    --trace_dir imputation_traces_1stage \
    -resume
