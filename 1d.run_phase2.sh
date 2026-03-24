#!/bin/bash
# Local test run: Phase 2 — SV imputation using Phase 1 SNP-imputed data
# Input: per-batch SNP-imputed VCFs from Phase 1 (imputation_output_phase1/snp_imputed_vcfs/)
set -euo pipefail

nextflow run main_phase2.nf \
    --phase1_snp_imputed_dir imputation_output_phase1/snp_imputed_vcfs \
    --sv_ref_panel /Users/thanhdng/Documents/databases/human/elife_2025_sv_panel/panel.888samples.full.vcf.gz \
    --plink_map_dir data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir imputation_output_phase2 \
    --trace_dir imputation_traces_phase2 \
    -resume
