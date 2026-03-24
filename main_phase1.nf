#!/usr/bin/env nextflow
/*
========================================================================================
                         sv-imputation — Phase 1: SNP Imputation
========================================================================================
  Inputs  : PLINK bfiles (GRCh37) + 1000 Genomes Project 3202-sample reference panel
  Outputs : Per-chromosome merged SNP-imputed VCFs (one file per chromosome,
            all input batches merged), ready for Phase 2 SV imputation
========================================================================================
*/

/*
 Define the default parameters
*/

// Input data paths (combined directory + bfiles prefix)
params.input_manifest     = "$baseDir/data/input_manifest.csv"

// Reference genome FASTA files from Ensembl
params.genome_build_37    = "$baseDir/data/genomes/GRCh37.fa"
params.genome_build_38    = "$baseDir/data/genomes/GRCh38.fa"

// Chain file: https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
params.chain_file         = "$baseDir/data/genomes/GRCh37_to_GRCh38.chain.gz"

// 1000 Genomes Project 3202-sample SNP reference panel
// Expects per-chromosome files: chr{N}.vcf.gz (+ .csi index)
params.snp_ref_panel_dir  = "/Users/thanhdng/Documents/databases/human/HC_1000G_hg38/rsid_vcfs"

// PLINK-format genetic map directory
// Files: plink.chr{N}.GRCh38.map  (numeric chrom IDs, no "chr" prefix)
params.plink_map_dir      = "$baseDir/data/plink.GRCh38.map/no_chr_in_chrom_field"

// Other parameters
params.chromosomes        = 1..22
params.outdir             = "./imputation_output_phase1"
params.trace_dir          = "./imputation_traces_phase1"


log.info """\
================================================================
              nf-imputation-sv  —  Phase 1: SNP Imputation
================================================================
  Input manifest : ${params.input_manifest}
  SNP ref panel  : ${params.snp_ref_panel_dir}
  Chromosomes    : ${params.chromosomes}
  Output dir     : ${params.outdir}
================================================================
"""

nextflow.enable.dsl=2

workflow {

    ref37 = file(params.genome_build_37)
    ref38 = file(params.genome_build_38)
    chain = file(params.chain_file)

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

    // Partition by genome build
    g37_bfile_ch = bfile_ch.filter { base, build, bed, bim, fam -> build == 'g37' }
                           .map    { base, build, bed, bim, fam -> tuple(base, bed, bim, fam) }
    g38_bfile_ch = bfile_ch.filter { base, build, bed, bim, fam -> build == 'g38' }
                           .map    { base, build, bed, bim, fam -> tuple(base, bed, bim, fam) }

    // Step 1+2: g37 batches → convert (GRCh37 ref) → liftover to GRCh38
    g37_vcf_ch    = ArrayBfile_to_VCF(g37_bfile_ch, ref37)
    g37_lifted_ch = LiftoverVcf(g37_vcf_ch, chain, ref37, ref38)
        .map { batch, vcf, csi, rejected -> tuple(batch, vcf, csi) }

    // Step 1 only: g38 batches → convert directly (GRCh38 ref, no liftover needed)
    g38_lifted_ch = ArrayBfile_to_VCF_38(g38_bfile_ch, ref38)

    // Merge both streams — from here on all batches are GRCh38
    lifted_ch = g37_lifted_ch.mix(g38_lifted_ch)

    // Step 3: Collect union of array positions (CHROM/POS) across all batches
    array_pos_ch = CollectArrayPositions( lifted_ch.map { batch, vcf, csi -> vcf }.collect() )

    // Step 4: Split lifted array VCFs by chromosome
    split_array_ch = SplitArrayVcf( lifted_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )
    // emits: (chrom, batch, vcf, csi)

    // Per-chromosome genetic map files (staged as process inputs)
    map_files_ch = Channel.fromList(params.chromosomes.toList())
        .map { chrom -> tuple(chrom, file("${params.plink_map_dir}/plink.chr${chrom}.GRCh38.map")) }

    // Raw 1000G per-chromosome VCFs
    snp_raw_ref_ch = Channel.fromList(params.chromosomes.toList())
        .map { chrom ->
            tuple(chrom,
                  file("${params.snp_ref_panel_dir}/chr${chrom}.vcf.gz"),
                  file("${params.snp_ref_panel_dir}/chr${chrom}.vcf.gz.csi"))
        }
    // emits: (chrom, raw_vcf, raw_csi)

    // Step 5a: Filter 1000G panel to MAF > 5% SNPs only (per chromosome)
    maf_ref_ch = FilterRefByMAF( snp_raw_ref_ch )
    // emits: (chrom, maf_vcf, maf_csi)

    // Step 5b: Filter 1000G panel to the union of array variant positions (per chromosome)
    array_ref_ch = FilterRefByArrayPositions( array_pos_ch.combine(snp_raw_ref_ch) )
    // emits: (chrom, array_vcf, array_csi)

    // Step 5c: Merge both filtered sets and remove duplicates → final SNP ref panel
    snp_ref_ch = MergeAndDedup( maf_ref_ch.join(array_ref_ch, by: 0) )
    // emits: (chrom, snp_ref_vcf, snp_ref_csi)

    // Step 6: Impute common SNPs per batch × chromosome
    snp_imputed_ch = BeagleSNPImpute( split_array_ch.combine(snp_ref_ch, by: 0).combine(map_files_ch, by: 0) )
    // emits: (chrom, batch, vcf, csi)

    // Step 7: Merge all batches per chromosome
    snp_merge_input_ch = snp_imputed_ch
        .map { chrom, batch, vcf, csi -> tuple(chrom, vcf, csi) }
        .groupTuple(by: 0)
    MergeSNPBatches(snp_merge_input_ch)
    // emits: (chrom, merged_vcf, merged_csi)  — written to outdir/snp_imputed_merged/
}

// ─────────────────────────────────────────────────────────────────────────────
//  Reference panel construction (Steps 5a / 5b / 5c)
// ─────────────────────────────────────────────────────────────────────────────

process FilterRefByMAF {

    tag "chr${chrom}"
    publishDir "${params.outdir}/snp_ref_maf", mode: 'link'

    input:
    tuple val(chrom), path(snp_vcf), path(snp_csi)

    output:
    tuple val(chrom), path("chr${chrom}.maf5.vcf.gz"), path("chr${chrom}.maf5.vcf.gz.csi")

    cpus 2
    memory '8 GB'

    script:
    """
    # Keep only SNPs with MAF > 5% across all 1000G samples
    bcftools view \\
        --type snps \\
        --min-af 0.05:minor \\
        --output-type z \\
        --output chr${chrom}.maf5.vcf.gz \\
        --threads ${task.cpus} \\
        ${snp_vcf}

    bcftools index -f chr${chrom}.maf5.vcf.gz
    """
}

process FilterRefByArrayPositions {

    tag "chr${chrom}"
    publishDir "${params.outdir}/snp_ref_array", mode: 'link'

    input:
    tuple path(array_positions), val(chrom), path(snp_vcf), path(snp_csi)

    output:
    tuple val(chrom), path("chr${chrom}.array_variants.vcf.gz"), path("chr${chrom}.array_variants.vcf.gz.csi")

    cpus 2
    memory '8 GB'

    script:
    """
    # Build a targets file for this chromosome from the union of array positions
    awk -v c="${chrom}" '\$1 == c {print \$1"\t"\$2}' ${array_positions} > array_targets.txt

    # Extract those positions from the 1000G panel (any variant type, regardless of MAF)
    bcftools view \\
        --targets-file array_targets.txt \\
        --output-type z \\
        --output chr${chrom}.array_variants.vcf.gz \\
        --threads ${task.cpus} \\
        ${snp_vcf}

    bcftools index -f chr${chrom}.array_variants.vcf.gz

    rm -f array_targets.txt
    """
}

process MergeAndDedup {

    tag "chr${chrom}"
    publishDir "${params.outdir}/snp_ref_by_chr", mode: 'link'

    input:
    tuple val(chrom),
          path(maf_vcf),   path(maf_csi),
          path(array_vcf), path(array_csi)

    output:
    tuple val(chrom), path("chr${chrom}.snp_ref.vcf.gz"), path("chr${chrom}.snp_ref.vcf.gz.csi")

    cpus 2
    memory '8 GB'

    script:
    """
    # Union of MAF>5% SNPs and array-position variants, sorted and deduplicated
    bcftools concat \\
        --allow-overlaps \\
        --threads ${task.cpus} \\
        ${maf_vcf} ${array_vcf} \\
    | bcftools sort \\
    | bcftools norm \\
        --rm-dup all \\
        -Oz -o chr${chrom}.snp_ref.vcf.gz

    bcftools index -f chr${chrom}.snp_ref.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  Imputation and batch merging
// ─────────────────────────────────────────────────────────────────────────────

process BeagleSNPImpute {

    tag "${batch}_chr${chrom}"
    publishDir "${params.outdir}/snp_imputed_vcfs", mode: 'copy'

    input:
    tuple val(chrom), val(batch), path(vcf), path(csi),
          path(snp_ref_vcf), path(snp_ref_csi),
          path(map_file)

    output:
    tuple val(chrom), val(batch),
          path("${batch}.chr${chrom}.snp_imputed.vcf.gz"),
          path("${batch}.chr${chrom}.snp_imputed.vcf.gz.csi")

    cpus 4
    memory { 32.GB * Math.pow(2, task.attempt - 1) }
    errorStrategy { task.exitStatus in [137, 140, 143, 247] ? 'retry' : 'finish' }
    maxRetries 3

    script:
    def mem_gb = (task.memory.toGiga() - 1) as int
    """
    java -Xmx${mem_gb}g -jar ${projectDir}/bin/beagle.27Feb25.75f.jar \\
        gt=${vcf} \\
        ref=${snp_ref_vcf} \\
        map=${map_file} \\
        chrom=${chrom} \\
        out=${batch}.chr${chrom}.snp_imputed \\
        nthreads=${task.cpus} \\
        ap=true \\
        gp=true

    bcftools index -f ${batch}.chr${chrom}.snp_imputed.vcf.gz
    """
}

process MergeSNPBatches {

    tag "chr${chrom}"
    publishDir "${params.outdir}/snp_imputed_merged", mode: 'copy'

    input:
    tuple val(chrom), path(vcfs), path(csis)

    output:
    tuple val(chrom), path("chr${chrom}.snp_merged.vcf.gz"), path("chr${chrom}.snp_merged.vcf.gz.csi")

    cpus 4
    memory '8 GB'

    script:
    def vcf_list = (vcfs instanceof List ? vcfs : [vcfs]).join(' ')
    def n        = (vcfs instanceof List ? vcfs : [vcfs]).size()
    """
    if [ ${n} -eq 1 ]; then
        cp ${vcf_list} chr${chrom}.snp_merged.vcf.gz
    else
        bcftools merge \\
            --no-version \\
            --output-type z \\
            --output chr${chrom}.snp_merged.vcf.gz \\
            --threads ${task.cpus} \\
            ${vcf_list}
    fi

    bcftools index -f chr${chrom}.snp_merged.vcf.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  Array data preprocessing (Steps 1 – 4)
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
    # Step 1: Merge PAR into X (requires --make-bed, cannot combine with --recode)
    plink \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --merge-x no-fail \\
        --keep-allele-order \\
        --make-bed \\
        --memory 4000 \\
        --out ${batch}.mergedx

    # Step 2: Recode merged bfiles to VCF
    plink \\
        --bed ${batch}.mergedx.bed \\
        --bim ${batch}.mergedx.bim \\
        --fam ${batch}.mergedx.fam \\
        --keep-allele-order \\
        --recode vcf-iid \\
        --memory 4000 \\
        --out ${batch}

    # Rename PLINK numeric codes to Ensembl notation (23→X, 24→Y, 26→MT)
    printf '23\tX\n24\tY\n26\tMT\n' > chr_rename.txt

    bgzip -c ${batch}.vcf \\
    | bcftools annotate \\
        --rename-chrs chr_rename.txt \\
        -Oz -o ${batch}.renamed.vcf.gz
    bcftools index -f ${batch}.renamed.vcf.gz

    rm -f ${batch}.vcf ${batch}.mergedx.bed ${batch}.mergedx.bim \\
          ${batch}.mergedx.fam ${batch}.mergedx.nosex ${batch}.mergedx.log

    # Normalize: fix REF/ALT, remove duplicates
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
    path ref   // GRCh38 — used directly for normalisation (no liftover needed)

    output:
    tuple val(batch), path("${batch}.norm.vcf.gz"), path("${batch}.norm.vcf.gz.csi")

    cpus 2
    memory '4 GB'

    script:
    """
    # Step 1: Merge PAR into X
    plink \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --merge-x no-fail \\
        --keep-allele-order \\
        --make-bed \\
        --memory 4000 \\
        --out ${batch}.mergedx

    # Step 2: Recode to VCF
    plink \\
        --bed ${batch}.mergedx.bed \\
        --bim ${batch}.mergedx.bim \\
        --fam ${batch}.mergedx.fam \\
        --keep-allele-order \\
        --recode vcf-iid \\
        --memory 4000 \\
        --out ${batch}

    # Rename PLINK numeric codes to Ensembl notation
    printf '23\tX\n24\tY\n26\tMT\n' > chr_rename.txt

    bgzip -c ${batch}.vcf \\
    | bcftools annotate \\
        --rename-chrs chr_rename.txt \\
        -Oz -o ${batch}.renamed.vcf.gz
    bcftools index -f ${batch}.renamed.vcf.gz

    rm -f ${batch}.vcf ${batch}.mergedx.bed ${batch}.mergedx.bim \\
          ${batch}.mergedx.fam ${batch}.mergedx.nosex ${batch}.mergedx.log

    # Normalize against GRCh38 directly (input is already GRCh38)
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
