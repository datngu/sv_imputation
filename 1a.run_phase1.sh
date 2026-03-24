#!/bin/bash
# Local test run: Phase 1 — SNP imputation of array batches
# Output: per-batch × per-chromosome SNP-imputed VCFs ready for Phase 2
set -euo pipefail

nextflow run main_phase1.nf \
    --input_manifest data/input_manifest.csv \
    --genome_build_37 data/genomes/GRCh37.fa \
    --genome_build_38 data/genomes/GRCh38.fa \
    --chain_file data/genomes/GRCh37_to_GRCh38.chain.gz \
    --snp_ref_panel_dir /Users/thanhdng/Documents/databases/human/HC_1000G_hg38/rsid_vcfs \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir imputation_output_phase1 \
    --trace_dir imputation_traces_phase1 \
    -resume
