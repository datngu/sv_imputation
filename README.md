# sv-imputation

A Nextflow DSL2 pipeline for imputing structural variants (SVs) into SNP array data using Beagle 5.

The core idea is to use a long-read sequencing-derived SV reference panel and impute SV genotypes into
short-read / SNP array samples. The pipeline follows the approach described in
[elifesciences.org/reviewed-preprints/106115](https://elifesciences.org/reviewed-preprints/106115).

---

## Dependencies

| Tool | Version tested | Install |
|------|---------------|--------|
| [Nextflow](https://www.nextflow.io/) | ≥ 24.x | `curl -s https://get.nextflow.io \| bash` |
| [bcftools](https://samtools.github.io/bcftools/) | 1.22 | `conda install -c bioconda bcftools` |
| [htslib](https://github.com/samtools/htslib) (bgzip, tabix) | 1.22 | bundled with bcftools |
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | v1.9.0-b.8 | `conda install -c bioconda plink` |
| [Python](https://www.python.org/) | 3.10+ | system / conda |
| Java | OpenJDK 21 | `conda install -c conda-forge openjdk` |
| Picard (bundled) | 3.4.0 | `bin/picard.jar` |
| Beagle 5 (bundled) | 27Feb25.75f | `bin/beagle.27Feb25.75f.jar` |

### Quick conda environment

```bash
conda create -n sv_imputation -c bioconda -c conda-forge \
    bcftools plink openjdk python
conda activate sv_imputation
```

---

## Required input data

| File | Description |
|------|-------------|
| `data/array_plink/{batch}.bed/bim/fam` | PLINK 1 bfiles, one prefix per array batch |
| `data/input_manifest.csv` | CSV listing bfile prefixes (see below) |
| `data/genomes/GRCh37.fa` + `.fai` | GRCh37 reference FASTA |
| `data/genomes/GRCh38.fa` + `.fai` | GRCh38 reference FASTA |
| `data/genomes/GRCh37_to_GRCh38.chain.gz` | Liftover chain file |
| `data/ref_panel/panel.888samples.full.vcf.gz` + `.csi` | Full SV reference panel (Ensembl chr convention) |

### Indexing reference files

```bash
# Reference FASTA
samtools faidx data/genomes/GRCh37.fa
samtools faidx data/genomes/GRCh38.fa

# Reference panel VCF (CSI)
bcftools index -f data/ref_panel/panel.888samples.full.vcf.gz
```

### Input manifest (`data/input_manifest.csv`)

One row per array batch. The `prefix` is the relative path to the PLINK bfile set
(without `.bed`/`.bim`/`.fam` extension):

```
prefix
data/array_plink/batch1
data/array_plink/batch2
```

---

## Workflow

```
Step 1  ArrayBfile_to_VCF      PLINK bfiles → normalised VCF (GRCh37, autosomes)
Step 2  LiftoverVcf            GRCh37 → GRCh38 (Picard LiftoverVcf)
Step 3  CollectArrayPositions  Pool CHROM/POS from all batches into one positions file
        FilterRefPanel         Filter full ref panel: keep array positions + all Sniffles2 SVs
Step 4  SplitRefPanel          Split reduced ref panel by chromosome
        SplitArrayVcf          Split each lifted batch VCF by chromosome
Step 5  BeagleImpute           Beagle 5 phasing + imputation, per batch × per chromosome
```

---

## Running the pipeline

```bash
bash 1.run_test.sh
```

Or manually:

```bash
nextflow run main.nf \
    --input_manifest    data/input_manifest.csv \
    --genome_build_37   data/genomes/GRCh37.fa \
    --genome_build_38   data/genomes/GRCh38.fa \
    --chain_file        data/genomes/GRCh37_to_GRCh38.chain.gz \
    --impute_ref_panel  data/ref_panel/panel.888samples.full.vcf.gz \
    --chromosomes       1-22 \
    --outdir            imputation_output \
    --trace_dir         imputation_traces \
    -resume
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_manifest` | `data/input_manifest.csv` | CSV of bfile prefixes |
| `--genome_build_37` | `data/genomes/GRCh37.fa` | GRCh37 FASTA |
| `--genome_build_38` | `data/genomes/GRCh38.fa` | GRCh38 FASTA |
| `--chain_file` | `data/genomes/GRCh37_to_GRCh38.chain.gz` | Liftover chain |
| `--impute_ref_panel` | `data/ref_panel/panel.888samples.full.vcf.gz` | Full reference panel VCF |
| `--chromosomes` | `1-22` | Chromosomes to process (Ensembl convention) |
| `--outdir` | `./imputation_output` | Output directory |
| `--trace_dir` | `./imputation_traces` | Nextflow trace / logs directory |

---

## Outputs

```
imputation_output/
    array_vcfs/               {batch}.norm.vcf.gz + .csi      (GRCh37, normalised)
    lifted_vcfs/              {batch}.ref38.vcf.gz + .csi     (GRCh38 lifted)
                              {batch}.rejected.vcf.gz          (liftover rejects)
    reduced_ref_panel/        reduced_ref.vcf.gz + .csi        (SVs + array positions)
    ref_panel_by_chr/         chr{N}.ref.vcf.gz + .csi         (per-chromosome ref)
    imputed_vcfs/             {batch}.chr{N}.imputed.vcf.gz + .csi
```

> All VCF indices use **CSI format** (`.csi`), produced by `bcftools index`.


## Alexey - meeting notes

- 2 folds imputations
- 1st fold: impute SVs into array data (using larger reference panel with many samples, but ony overlap array positions + SNPs with higher MAF (AC > 10), can be create a customized ref-panel or exlcude variants from the full panel during imputation)

- 2nd fold: impute SVs into array data (using smaller reference panel with fewer samples, but all SVs)