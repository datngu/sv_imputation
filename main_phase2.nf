#!/usr/bin/env nextflow
/*
========================================================================================
                         sv-imputation — Phase 2: SV Imputation
========================================================================================
  Inputs  : Per-batch SNP-imputed VCFs from Phase 1
            (snp_imputed_vcfs/{batch}.chr{N}.snp_imputed.vcf.gz)
          + Full SV reference panel (e.g. panel.888samples.full.vcf.gz)
  Process : Build a reduced ref panel containing the Phase-1-imputed SNP positions
            plus all SVs from the full ref → per-batch Beagle SV imputation → merge
  Outputs : Per-chromosome merged SV-imputed VCFs, one file per chromosome with
            all input batches merged
========================================================================================
*/

/*
 Define the default parameters
*/

// Directory containing Phase 1 per-batch SNP-imputed VCFs
// Expected file naming: {batch}.chr{N}.snp_imputed.vcf.gz (+ .csi)
params.phase1_snp_imputed_dir = "./imputation_output_phase1/snp_imputed_vcfs"

// Full SV reference panel VCF (all chromosomes, all samples)
params.sv_ref_panel           = "$baseDir/data/ref_panel/panel.888samples.full.vcf.gz"

// PLINK-format genetic map directory
// Files: plink.chr{N}.GRCh38.map  (numeric chrom IDs, no "chr" prefix)
params.plink_map_dir          = "$baseDir/data/plink.GRCh38.map/no_chr_in_chrom_field"

// Other parameters
params.chromosomes            = 1..22
params.outdir                 = "./imputation_output_phase2"
params.trace_dir              = "./imputation_traces_phase2"


log.info """\
================================================================
              nf-imputation-sv  —  Phase 2: SV Imputation
================================================================
  Phase 1 SNP-imputed dir : ${params.phase1_snp_imputed_dir}
  SV ref panel            : ${params.sv_ref_panel}
  Chromosomes             : ${params.chromosomes}
  Output dir              : ${params.outdir}
================================================================
"""

nextflow.enable.dsl=2

workflow {

    sv_ref_panel = file(params.sv_ref_panel)

    // Per-chromosome genetic map files (staged as process inputs)
    map_files_ch = Channel.fromList(params.chromosomes.toList())
        .map { chrom -> tuple(chrom, file("${params.plink_map_dir}/plink.chr${chrom}.GRCh38.map")) }

    // Read Phase 1 per-batch SNP-imputed VCFs
    // Filename pattern: {batch}.chr{N}.snp_imputed.vcf.gz
    snp_imputed_ch = Channel
        .fromPath("${params.phase1_snp_imputed_dir}/*.snp_imputed.vcf.gz")
        .map { vcf ->
            def m = (vcf.getName() =~ /^(.+)\.chr(\d+)\.snp_imputed\.vcf\.gz$/)
            if (!m) error "Unexpected Phase 1 VCF filename: ${vcf.getName()}"
            tuple(m[0][2].toInteger(), m[0][1], vcf, file("${vcf}.csi"))
        }
        .filter { chrom, batch, vcf, csi -> chrom in params.chromosomes.toList() }
    // emits: (chrom, batch, vcf, csi)

    // Step 1: Collect SNP positions from Phase 1 imputed VCFs.
    //         All batches share the same reference panel → identical positions per chromosome.
    //         Take one VCF per chromosome (first batch) to avoid redundant work.
    snp_pos_ch = CollectSNPPositions(
        snp_imputed_ch
            .groupTuple(by: 0)
            .map { chrom, batches, vcfs, csis -> vcfs instanceof List ? vcfs[0] : vcfs }
            .collect()
    )
    // emits: path("all_snp_positions.tsv")

    // Step 2: Build reduced SV reference panel:
    //         keep variants at Phase-1-imputed SNP positions + all SVs (Sniffles2 IDs)
    reduced_ref_ch = FilterRefPanel(snp_pos_ch, sv_ref_panel)
    // emits: (reduced_ref.vcf.gz, reduced_ref.vcf.gz.csi)

    // Step 3: Split reduced ref panel by chromosome
    split_ref_ch = SplitRefPanel( reduced_ref_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )
    // emits: (chrom, chr{N}.ref.vcf.gz, csi)

    // Step 4: SV imputation — per batch × chromosome (keep batches separate to manage memory)
    sv_imputed_ch = BeagleSVImpute(
        snp_imputed_ch.combine(split_ref_ch, by: 0).combine(map_files_ch, by: 0)
    )
    // emits: (chrom, batch, sv_imputed.vcf.gz, csi)

    // Step 5: Merge all batches per chromosome → final SV-imputed VCF
    sv_merge_input_ch = sv_imputed_ch
        .map { chrom, batch, vcf, csi -> tuple(chrom, vcf, csi) }
        .groupTuple(by: 0)
    MergeSVBatches(sv_merge_input_ch)
    // emits: (chrom, sv_merged.vcf.gz, csi) — written to outdir/sv_imputed_merged/
}

// ─────────────────────────────────────────────────────────────────────────────
//  Reference panel construction (Steps 1 – 3)
// ─────────────────────────────────────────────────────────────────────────────

process CollectSNPPositions {

    publishDir "${params.trace_dir}/snp_positions", mode: 'link'

    input:
    path vcfs   // all Phase 1 per-batch imputed VCFs collected into one list

    output:
    path "all_snp_positions.tsv"

    cpus 1
    memory '4 GB'

    script:
    """
    for vcf in ${vcfs}; do
        bcftools query -f '%CHROM\t%POS\n' "\$vcf"
    done | sort -k1,1 -k2,2n | uniq > all_snp_positions.tsv
    """
}

process FilterRefPanel {

    publishDir "${params.outdir}/reduced_ref_panel", mode: 'link'

    input:
    path positions
    path ref_panel

    output:
    tuple path("reduced_ref.vcf.gz"), path("reduced_ref.vcf.gz.csi")

    cpus 1
    memory '8 GB'

    script:
    """
    python3 ${projectDir}/bin/filter_ref_panel.py \\
        --positions ${positions} \\
        --ref ${ref_panel} \\
        | bgzip -c > reduced_ref.vcf.gz

    bcftools index -f reduced_ref.vcf.gz
    """
}

process SplitRefPanel {

    tag "chr${chrom}"
    publishDir "${params.outdir}/ref_panel_by_chr", mode: 'link'

    input:
    tuple path(ref_vcf), path(ref_csi), val(chrom)

    output:
    tuple val(chrom), path("chr${chrom}.ref.vcf.gz"), path("chr${chrom}.ref.vcf.gz.csi")

    cpus 2
    memory '4 GB'

    script:
    """
    bcftools view \\
        --regions ${chrom} \\
        --output-type z \\
        --output chr${chrom}.ref.vcf.gz \\
        ${ref_vcf}

    bcftools index -f chr${chrom}.ref.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  SV imputation and batch merging (Steps 4 – 5)
// ─────────────────────────────────────────────────────────────────────────────

process BeagleSVImpute {

    tag "${batch}_chr${chrom}"
    publishDir "${params.outdir}/sv_imputed_vcfs", mode: 'copy'

    input:
    tuple val(chrom), val(batch), path(vcf), path(csi),
          path(sv_ref_vcf), path(sv_ref_csi),
          path(map_file)

    output:
    tuple val(chrom), val(batch),
          path("${batch}.chr${chrom}.sv_imputed.vcf.gz"),
          path("${batch}.chr${chrom}.sv_imputed.vcf.gz.csi")

    cpus 4
    memory { 32.GB * Math.pow(2, task.attempt - 1) }
    errorStrategy { task.exitStatus in [137, 140, 143, 247] ? 'retry' : 'finish' }
    maxRetries 3

    script:
    def mem_gb = (task.memory.toGiga() - 1) as int
    """
    java -Xmx${mem_gb}g -jar ${projectDir}/bin/beagle.27Feb25.75f.jar \\
        gt=${vcf} \\
        ref=${sv_ref_vcf} \\
        map=${map_file} \\
        chrom=${chrom} \\
        out=${batch}.chr${chrom}.sv_imputed \\
        nthreads=${task.cpus} \\
        ap=true \\
        gp=true

    bcftools index -f ${batch}.chr${chrom}.sv_imputed.vcf.gz
    """
}

process MergeSVBatches {

    tag "chr${chrom}"
    publishDir "${params.outdir}/sv_imputed_merged", mode: 'copy'

    input:
    tuple val(chrom), path(vcfs), path(csis)

    output:
    tuple val(chrom), path("chr${chrom}.sv_merged.vcf.gz"), path("chr${chrom}.sv_merged.vcf.gz.csi")

    cpus 4
    memory '8 GB'

    script:
    def vcf_list = (vcfs instanceof List ? vcfs : [vcfs]).join(' ')
    def n        = (vcfs instanceof List ? vcfs : [vcfs]).size()
    """
    if [ ${n} -eq 1 ]; then
        cp ${vcf_list} chr${chrom}.sv_merged.vcf.gz
    else
        bcftools merge \\
            --no-version \\
            --output-type z \\
            --output chr${chrom}.sv_merged.vcf.gz \\
            --threads ${task.cpus} \\
            ${vcf_list}
    fi

    bcftools index -f chr${chrom}.sv_merged.vcf.gz
    """
}
