#!/bin/bash
# Test run: 1-stage imputation (SV imputation directly on array data, no SNP pre-imputation)
set -euo pipefail

nextflow run main.nf \
    --input_manifest data/input_manifest.csv \
    --genome_build_37 data/genomes/GRCh37.fa \
    --genome_build_38 data/genomes/GRCh38.fa \
    --chain_file data/genomes/GRCh37_to_GRCh38.chain.gz \
    --impute_ref_panel /Users/thanhdng/Documents/databases/human/elife_2025_sv_panel/panel.888samples.full.vcf.gz \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --imputation_mode 1stage \
    --outdir imputation_output_1stage \
    --trace_dir imputation_traces_1stage \
    -resume
