#!/usr/bin/env nextflow
/*
========================================================================================
                         sv-imputation — 1-stage: direct SV imputation
========================================================================================
  Inputs  : PLINK bfiles (GRCh37 or GRCh38) + full SV reference panel
  Process : Preprocess bfiles → (liftover if GRCh37) → collect array positions
            → reduce SV ref panel (array positions + all SVs) → split by chrom
            → per-batch Beagle SV imputation → merge batches per chromosome.
  Note    : No intermediate SNP imputation. Faster but lower SV imputation
            accuracy compared to the 2-stage approach (main_phase1.nf + main_phase2.nf).
========================================================================================
*/

/*
 Define the default parameters
*/

params.input_manifest  = "$baseDir/data/input_manifest.csv"

params.genome_build_37 = "$baseDir/data/genomes/GRCh37.fa"
params.genome_build_38 = "$baseDir/data/genomes/GRCh38.fa"

// Chain file: https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
params.chain_file      = "$baseDir/data/genomes/GRCh37_to_GRCh38.chain.gz"

// Full SV reference panel VCF (all chromosomes, all samples)
params.sv_ref_panel    = "$baseDir/data/ref_panel/panel.888samples.full.vcf.gz"

// PLINK-format genetic map directory
// Files: plink.chr{N}.GRCh38.map  (numeric chrom IDs, no "chr" prefix)
params.plink_map_dir   = "$baseDir/data/plink.GRCh38.map/no_chr_in_chrom_field"

params.chromosomes     = 1..22
params.outdir          = "./imputation_output_1stage"
params.trace_dir       = "./imputation_traces_1stage"


log.info """\
================================================================
              nf-imputation-sv  —  1-stage SV imputation
================================================================
  Input manifest : ${params.input_manifest}
  SV ref panel   : ${params.sv_ref_panel}
  Chromosomes    : ${params.chromosomes}
  Output dir     : ${params.outdir}
================================================================
"""

nextflow.enable.dsl=2

workflow {

    ref37      = file(params.genome_build_37)
    ref38      = file(params.genome_build_38)
    chain      = file(params.chain_file)
    sv_ref     = file(params.sv_ref_panel)

    // Per-chromosome genetic map files (staged as process inputs)
    map_files_ch = Channel.fromList(params.chromosomes.toList())
        .map { chrom -> tuple(chrom, file("${params.plink_map_dir}/plink.chr${chrom}.GRCh38.map")) }

    // Read manifest — partition by genome build
    bfile_ch = Channel
        .fromPath(params.input_manifest)
        .splitCsv(header: true)
        .filter { row -> row.prefix }
        .map { row ->
            def prefix = row.prefix
            def base   = prefix.tokenize('/')[-1]
            def build  = row.genome_build ?: 'g37'
            tuple(base, build, file("${prefix}.bed"), file("${prefix}.bim"), file("${prefix}.fam"))
        }

    g37_bfile_ch = bfile_ch.filter { base, build, bed, bim, fam -> build == 'g37' }
                           .map    { base, build, bed, bim, fam -> tuple(base, bed, bim, fam) }
    g38_bfile_ch = bfile_ch.filter { base, build, bed, bim, fam -> build == 'g38' }
                           .map    { base, build, bed, bim, fam -> tuple(base, bed, bim, fam) }

    // Step 1: Convert bfiles to normalised VCF
    // Step 2: Liftover g37 batches → GRCh38; g38 batches pass through directly
    g37_vcf_ch    = ArrayBfile_to_VCF(g37_bfile_ch, ref37)
    g37_lifted_ch = LiftoverVcf(g37_vcf_ch, chain, ref37, ref38)
        .map { batch, vcf, csi, rejected -> tuple(batch, vcf, csi) }

    g38_lifted_ch = ArrayBfile_to_VCF_38(g38_bfile_ch, ref38)

    // All batches are now GRCh38
    lifted_ch = g37_lifted_ch.mix(g38_lifted_ch)

    // Step 3: Collect union of array positions (CHROM/POS) across all batches
    array_pos_ch = CollectArrayPositions( lifted_ch.map { batch, vcf, csi -> vcf }.collect() )

    // Step 4: Reduce SV ref panel: keep variants at array positions + all SVs (Sniffles2 IDs)
    reduced_ref_ch = FilterRefPanel(array_pos_ch, sv_ref)
    // emits: (reduced_ref.vcf.gz, csi)

    // Step 5: Split reduced SV ref panel by chromosome
    split_ref_ch = SplitRefPanel( reduced_ref_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )
    // emits: (chrom, chr{N}.ref.vcf.gz, csi)

    // Step 6: Split lifted array VCFs by chromosome (per batch)
    split_array_ch = SplitArrayVcf( lifted_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )
    // emits: (chrom, batch, vcf, csi)

    // Step 7: SV imputation — per batch × chromosome
    sv_imputed_ch = BeagleSVImpute(
        split_array_ch.combine(split_ref_ch, by: 0).combine(map_files_ch, by: 0)
    )
    // emits: (chrom, batch, sv_imputed.vcf.gz, csi)

    // Step 8: Merge all batches per chromosome
    sv_merge_input_ch = sv_imputed_ch
        .map { chrom, batch, vcf, csi -> tuple(chrom, vcf, csi) }
        .groupTuple(by: 0)
    MergeSVBatches(sv_merge_input_ch)
    // emits: (chrom, sv_merged.vcf.gz, csi) — written to outdir/sv_imputed_merged/
}

// ─────────────────────────────────────────────────────────────────────────────
//  Array data preprocessing (Steps 1 – 3)
// ─────────────────────────────────────────────────────────────────────────────

process ArrayBfile_to_VCF {

    tag "${batch}"
    publishDir "${params.outdir}/array_vcfs", mode: 'copy'

    input:
    tuple val(batch), path(bed), path(bim), path(fam)
    path ref   // GRCh37 — used for normalisation before liftover

    output:
    tuple val(batch), path("${batch}.norm.vcf.gz"), path("${batch}.norm.vcf.gz.csi")

    cpus 2
    memory '4 GB'

    script:
    """
    plink \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --merge-x no-fail \\
        --keep-allele-order \\
        --make-bed \\
        --memory 4000 \\
        --out ${batch}.mergedx

    plink \\
        --bed ${batch}.mergedx.bed \\
        --bim ${batch}.mergedx.bim \\
        --fam ${batch}.mergedx.fam \\
        --keep-allele-order \\
        --recode vcf-iid \\
        --memory 4000 \\
        --out ${batch}

    printf '23\tX\n24\tY\n26\tMT\n' > chr_rename.txt

    bgzip -c ${batch}.vcf \\
    | bcftools annotate \\
        --rename-chrs chr_rename.txt \\
        -Oz -o ${batch}.renamed.vcf.gz
    bcftools index -f ${batch}.renamed.vcf.gz

    rm -f ${batch}.vcf ${batch}.mergedx.bed ${batch}.mergedx.bim \\
          ${batch}.mergedx.fam ${batch}.mergedx.nosex ${batch}.mergedx.log

    bcftools norm \\
        -c sx \\
        --rm-dup all \\
        -f ${ref} \\
        -Oz \\
        -o ${batch}.norm.vcf.gz \\
        ${batch}.renamed.vcf.gz

    rm -f ${batch}.renamed.vcf.gz ${batch}.renamed.vcf.gz.csi

    bcftools index -f ${batch}.norm.vcf.gz

    bcftools norm -c e -f ${ref} ${batch}.norm.vcf.gz 2>&1 | grep -i "mismatch" || true
    """
}

process ArrayBfile_to_VCF_38 {

    tag "${batch}"
    publishDir "${params.outdir}/array_vcfs", mode: 'copy'

    input:
    tuple val(batch), path(bed), path(bim), path(fam)
    path ref   // GRCh38 — normalise directly, no liftover needed

    output:
    tuple val(batch), path("${batch}.norm.vcf.gz"), path("${batch}.norm.vcf.gz.csi")

    cpus 2
    memory '4 GB'

    script:
    """
    plink \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --merge-x no-fail \\
        --keep-allele-order \\
        --make-bed \\
        --memory 4000 \\
        --out ${batch}.mergedx

    plink \\
        --bed ${batch}.mergedx.bed \\
        --bim ${batch}.mergedx.bim \\
        --fam ${batch}.mergedx.fam \\
        --keep-allele-order \\
        --recode vcf-iid \\
        --memory 4000 \\
        --out ${batch}

    printf '23\tX\n24\tY\n26\tMT\n' > chr_rename.txt

    bgzip -c ${batch}.vcf \\
    | bcftools annotate \\
        --rename-chrs chr_rename.txt \\
        -Oz -o ${batch}.renamed.vcf.gz
    bcftools index -f ${batch}.renamed.vcf.gz

    rm -f ${batch}.vcf ${batch}.mergedx.bed ${batch}.mergedx.bim \\
          ${batch}.mergedx.fam ${batch}.mergedx.nosex ${batch}.mergedx.log

    bcftools norm \\
        -c sx \\
        --rm-dup all \\
        -f ${ref} \\
        -Oz \\
        -o ${batch}.norm.vcf.gz \\
        ${batch}.renamed.vcf.gz

    rm -f ${batch}.renamed.vcf.gz ${batch}.renamed.vcf.gz.csi

    bcftools index -f ${batch}.norm.vcf.gz

    bcftools norm -c e -f ${ref} ${batch}.norm.vcf.gz 2>&1 | grep -i "mismatch" || true
    """
}

process LiftoverVcf {

    tag "${batch}"
    publishDir "${params.outdir}/lifted_vcfs", mode: 'copy'

    input:
    tuple val(batch), path(vcf), path(tbi)
    path chain
    path ref37
    path ref38

    output:
    tuple val(batch), path("${batch}.ref38.vcf.gz"), path("${batch}.ref38.vcf.gz.csi"), path("${batch}.rejected.vcf.gz")

    cpus 2
    memory '8 GB'

    script:
    """
    bcftools +liftover \\
        --no-version \\
        -Ou \\
        ${vcf} \\
        -- \\
        --src-fasta-ref ${ref37} \\
        --fasta-ref ${ref38} \\
        --chain ${chain} \\
        --reject ${batch}.rejected.vcf.gz \\
        --reject-type z \\
    | bcftools sort \\
        -Oz -o ${batch}.ref38.vcf.gz

    bcftools index -f ${batch}.ref38.vcf.gz
    """
}

process CollectArrayPositions {

    publishDir "${params.trace_dir}/array_positions", mode: 'link'

    input:
    path vcfs

    output:
    path "all_array_positions.tsv"

    cpus 1
    memory '4 GB'

    script:
    """
    for vcf in ${vcfs}; do
        bcftools query -f '%CHROM\t%POS\n' "\$vcf"
    done | sort -k1,1 -k2,2n | uniq > all_array_positions.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  SV reference panel construction (Steps 4 – 5)
// ─────────────────────────────────────────────────────────────────────────────

process FilterRefPanel {

    publishDir "${params.outdir}/reduced_ref_panel", mode: 'link'

    input:
    path positions
    path ref_panel

    output:
    tuple path("reduced_ref.vcf.gz"), path("reduced_ref.vcf.gz.csi")

    cpus 1
    memory '32 GB'

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
//  SV imputation and batch merging (Steps 6 – 8)
// ─────────────────────────────────────────────────────────────────────────────

process SplitArrayVcf {

    tag "${batch}_chr${chrom}"

    input:
    tuple val(batch), path(vcf), path(tbi), val(chrom)

    output:
    tuple val(chrom), val(batch), path("${batch}.chr${chrom}.vcf.gz"), path("${batch}.chr${chrom}.vcf.gz.csi")

    cpus 1
    memory '4 GB'

    script:
    """
    bcftools view \\
        --regions ${chrom} \\
        --output-type z \\
        --output ${batch}.chr${chrom}.vcf.gz \\
        ${vcf}

    bcftools index -f ${batch}.chr${chrom}.vcf.gz
    """
}

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