# sv-imputation

A Nextflow DSL2 pipeline for imputing structural variants (SVs) into SNP array data using Beagle 5.

The core idea is to use a long-read sequencing-derived SV reference panel and impute SV genotypes into
short-read / SNP array samples. The pipeline follows the approach described in
[elifesciences.org/reviewed-preprints/106115](https://elifesciences.org/reviewed-preprints/106115).

The pipeline is split into two independent phases:

- **Phase 1** (`main_phase1.nf`): Preprocess array bfiles, (optionally) liftover to GRCh38, and impute
  common SNPs per batch using a large 1000 Genomes reference panel. Outputs one
  SNP-imputed VCF per batch × per chromosome.
- **Phase 2** (`main_phase2.nf`): Use the Phase 1 per-batch SNP-imputed VCFs as input, build a reduced
  SV reference panel (SNP positions from Phase 1 + all SVs), and run Beagle SV imputation
  per batch × per chromosome. Merge batches per chromosome at the end.

---

## Dependencies

| Tool | Version tested | Install |
|------|---------------|--------|
| [Nextflow](https://www.nextflow.io/) | ≥ 24.x | `curl -s https://get.nextflow.io \| bash` |
| [bcftools](https://samtools.github.io/bcftools/) + liftover plugin | 1.22 | `conda install -c bioconda bcftools bcftools-liftover-plugin` |
| [htslib](https://github.com/samtools/htslib) (bgzip) | 1.22 | bundled with bcftools |
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | v1.9.0-b.8 | `conda install -c bioconda plink` |
| [Python](https://www.python.org/) | 3.12 | system / conda |
| Java | OpenJDK 17 | `conda install -c conda-forge openjdk=17` |
| Beagle 5 (bundled) | 27Feb25.75f | `bin/beagle.27Feb25.75f.jar` |

### Docker / Singularity image

A pre-built image with all tools is available on Docker Hub:

```
ndatth/sv-imputation:latest
```

Use the `sv_impute` Nextflow profile (see `nextflow.config`) to pull it automatically.

### Quick conda environment

```bash
conda create -n sv_imputation -c bioconda -c conda-forge \
    bcftools bcftools-liftover-plugin plink openjdk=17 python=3.12
conda activate sv_imputation
```

---

## Required input data

### Phase 1

| File / directory | Description |
|-----------------|-------------|
| `data/array_plink/{batch}.bed/bim/fam` | PLINK 1 bfiles, one prefix per array batch |
| `data/input_manifest.csv` | CSV with `prefix` and `genome_build` columns (see below) |
| `data/genomes/GRCh37.fa` + `.fai` | GRCh37 reference FASTA (for g37 batches) |
| `data/genomes/GRCh38.fa` + `.fai` | GRCh38 reference FASTA |
| `data/genomes/GRCh37_to_GRCh38.chain.gz` | Liftover chain file (for g37 batches) |
| `{snp_ref_panel_dir}/chr{N}.vcf.gz` + `.csi` | 1000G SNP ref panel, one VCF per chromosome |
| `data/plink.GRCh38.map/no_chr_in_chrom_field/plink.chr{N}.GRCh38.map` | Per-chromosome genetic maps |

### Phase 2

| File / directory | Description |
|-----------------|-------------|
| `imputation_output_phase1/snp_imputed_vcfs/` | Phase 1 output directory (per-batch SNP-imputed VCFs) |
| `data/ref_panel/panel.888samples.full.vcf.gz` + `.csi` | Full SV reference panel |
| `data/plink.GRCh38.map/no_chr_in_chrom_field/plink.chr{N}.GRCh38.map` | Same genetic maps as Phase 1 |

### Indexing reference files

```bash
samtools faidx data/genomes/GRCh37.fa
samtools faidx data/genomes/GRCh38.fa
bcftools index -f data/ref_panel/panel.888samples.full.vcf.gz
```

### Input manifest (`data/input_manifest.csv`)

One row per array batch. `prefix` is the path to the PLINK bfile set (without extension).
`genome_build` is `g37` (needs liftover) or `g38` (already GRCh38, no liftover):

```
prefix,genome_build
data/array_plink/batch1,g37
data/array_plink/batch2,g38
```

---

## Workflow overview

### Phase 1 — SNP imputation (`main_phase1.nf`)

```
Step 1  ArrayBfile_to_VCF       PLINK bfiles → normalised VCF (GRCh37 ref, g37 batches)
        ArrayBfile_to_VCF_38    PLINK bfiles → normalised VCF (GRCh38 ref, g38 batches)
Step 2  LiftoverVcf             GRCh37 → GRCh38 (bcftools +liftover, g37 batches only)
Step 3  CollectArrayPositions   Union of CHROM/POS across all batches
Step 4  SplitArrayVcf           Split each lifted VCF by chromosome
Step 5a FilterRefByMAF          1000G panel → SNPs with MAF > 5%
Step 5b FilterRefByArrayPositions 1000G panel → variants at array positions (any MAF)
Step 5c MergeAndDedup           Union of 5a + 5b → final per-chromosome SNP ref panel
Step 6  BeagleSNPImpute         Beagle 5 SNP imputation, per batch × per chromosome
Step 7  MergeSNPBatches         Merge all batches per chromosome
```

### Phase 2 — SV imputation (`main_phase2.nf`)

```
Step 1  CollectSNPPositions     Union of CHROM/POS from all Phase 1 per-batch VCFs
Step 2  FilterRefPanel          Reduce full SV ref: keep Phase-1 SNP positions + all SVs
Step 3  SplitRefPanel           Split reduced SV ref by chromosome
Step 4  BeagleSVImpute          Beagle 5 SV imputation, per batch × per chromosome
Step 5  MergeSVBatches          Merge all batches per chromosome
```

---

## Run scripts

| Script | Target | Description |
|--------|--------|-------------|
| `1c.run_phase1.sh` | local | Phase 1 — local test run |
| `1d.run_phase2.sh` | local | Phase 2 — local test run |
| `2c.run_saga_phase1.sh` | SAGA HPC | Phase 1 — SBATCH job (Singularity) |
| `2d.run_saga_phase2.sh` | SAGA HPC | Phase 2 — SBATCH job (Singularity) |


### Phase 1 (local)

```bash
bash 1c.run_phase1.sh
```

Or manually:

```bash
nextflow run main_phase1.nf \
    --input_manifest    data/input_manifest.csv \
    --genome_build_37   data/genomes/GRCh37.fa \
    --genome_build_38   data/genomes/GRCh38.fa \
    --chain_file        data/genomes/GRCh37_to_GRCh38.chain.gz \
    --snp_ref_panel_dir /path/to/HC_1000G_hg38/rsid_vcfs \
    --plink_map_dir     data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir            imputation_output_phase1 \
    --trace_dir         imputation_traces_phase1 \
    -resume
```

### Phase 2 (local)

```bash
bash 1d.run_phase2.sh
```

Or manually:

```bash
nextflow run main_phase2.nf \
    --phase1_snp_imputed_dir imputation_output_phase1/snp_imputed_vcfs \
    --sv_ref_panel      data/ref_panel/panel.888samples.full.vcf.gz \
    --plink_map_dir     data/plink.GRCh38.map/no_chr_in_chrom_field \
    --outdir            imputation_output_phase2 \
    --trace_dir         imputation_traces_phase2 \
    -resume
```

### SAGA HPC

```bash
sbatch 2c.run_saga_phase1.sh   # Phase 1
sbatch 2d.run_saga_phase2.sh   # Phase 2 (run after Phase 1 completes)
```

---

## Parameters

### Phase 1 (`main_phase1.nf`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_manifest` | `data/input_manifest.csv` | CSV with `prefix` and `genome_build` columns |
| `--genome_build_37` | `data/genomes/GRCh37.fa` | GRCh37 FASTA (for g37 batches) |
| `--genome_build_38` | `data/genomes/GRCh38.fa` | GRCh38 FASTA |
| `--chain_file` | `data/genomes/GRCh37_to_GRCh38.chain.gz` | Liftover chain |
| `--snp_ref_panel_dir` | *(required)* | Directory of 1000G per-chromosome VCFs |
| `--plink_map_dir` | `data/plink.GRCh38.map/no_chr_in_chrom_field` | Genetic map directory |
| `--chromosomes` | `1..22` | Chromosomes to process |
| `--outdir` | `./imputation_output_phase1` | Output directory |
| `--trace_dir` | `./imputation_traces_phase1` | Trace / log directory |

### Phase 2 (`main_phase2.nf`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--phase1_snp_imputed_dir` | `./imputation_output_phase1/snp_imputed_vcfs` | Phase 1 per-batch VCF directory |
| `--sv_ref_panel` | `data/ref_panel/panel.888samples.full.vcf.gz` | Full SV reference panel VCF |
| `--plink_map_dir` | `data/plink.GRCh38.map/no_chr_in_chrom_field` | Genetic map directory |
| `--chromosomes` | `1..22` | Chromosomes to process |
| `--outdir` | `./imputation_output_phase2` | Output directory |
| `--trace_dir` | `./imputation_traces_phase2` | Trace / log directory |

---

## Outputs

### Phase 1

```
imputation_output_phase1/
    array_vcfs/               {batch}.norm.vcf.gz + .csi       (GRCh38, normalised)
    lifted_vcfs/              {batch}.ref38.vcf.gz + .csi      (GRCh38 lifted, g37 batches)
                              {batch}.rejected.vcf.gz           (liftover rejects)
    snp_ref_maf/              chr{N}.maf5.vcf.gz + .csi        (1000G SNPs, MAF > 5%)
    snp_ref_array/            chr{N}.array_variants.vcf.gz + .csi (1000G at array positions)
    snp_ref_by_chr/           chr{N}.snp_ref.vcf.gz + .csi     (final union SNP ref panel)
    snp_imputed_vcfs/         {batch}.chr{N}.snp_imputed.vcf.gz + .csi
    snp_imputed_merged/       chr{N}.snp_merged.vcf.gz + .csi  (all batches merged per chr)
```

### Phase 2

```
imputation_output_phase2/
    reduced_ref_panel/        reduced_ref.vcf.gz + .csi        (SNPs + SVs)
    ref_panel_by_chr/         chr{N}.ref.vcf.gz + .csi         (per-chromosome SV ref)
    sv_imputed_vcfs/          {batch}.chr{N}.sv_imputed.vcf.gz + .csi
    sv_imputed_merged/        chr{N}.sv_merged.vcf.gz + .csi   (all batches merged per chr)
```

> All VCF indices use **CSI format** (`.csi`), produced by `bcftools index`.