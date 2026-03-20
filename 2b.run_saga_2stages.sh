#!/bin/bash
#SBATCH --job-name=sv_imputation_2stages
#SBATCH --account=nn9114k
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=sv_imputation_2stages_%j.log

# Load Nextflow
module load Nextflow/24.04.2

# Set Singularity cache and tmp directories
mkdir -p /cluster/projects/nn9114k/datngu/singularity/tmp

export SINGULARITY_CACHEDIR="/cluster/projects/nn9114k/datngu/singularity"
export SINGULARITY_TMPDIR="/cluster/projects/nn9114k/datngu/singularity/tmp"
export NXF_SINGULARITY_CACHEDIR="/cluster/projects/nn9114k/datngu/singularity"

nextflow run main.nf \
    -profile saga,singularity,sv_impute \
    --input_manifest data/input_manifest.csv \
    --genome_build_37 data/genomes/GRCh37.fa \
    --genome_build_38 data/genomes/GRCh38.fa \
    --chain_file data/genomes/GRCh37_to_GRCh38.chain.gz \
    --impute_ref_panel data/ref_panel/panel.888samples.full.vcf.gz \
    --snp_ref_panel_dir /path/to/HC_1000G_hg38/rsid_vcfs \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --imputation_mode 2stage \
    --outdir imputation_output_2stage \
    --trace_dir imputation_traces_2stage \
    -resume
