#!/usr/bin/env nextflow
/*
========================================================================================
                                 sv-imputation
========================================================================================
                    a nextflow pipeline for SV imputation
----------------------------------------------------------------------------------------
*/


/*
 Define the default parameters
*/ 

// Input data paths (combined directory + bfiles prefix)
params.input_manifest          = "$baseDir/data/input_manifest.csv"


// Reference genome FASTA files from Ensembl (ref37 and ref38 builds)
params.genome_build_37        = "$baseDir/data/genomes/GRCh37.fa"
params.genome_build_38        = "$baseDir/data/genomes/GRCh38.fa"

// Chain file: https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
params.chain_file              = "$baseDir/data/genomes/GRCh37_to_GRCh38.chain.gz"
params.impute_ref_panel        = "$baseDir/data/ref_panel/panel.888samples.full.vcf.gz"
// PLINK-format genetic map directory (per-chromosome files: plink.chr{N}.GRCh38.map)
// Use the no_chr_in_chrom_field variant — chrom column is numeric (1,2,...22),
// matching plink's output and the chromosomes used internally in this pipeline.
params.plink_map_dir           = "$baseDir/data/plink.GRCh38.map/no_chr_in_chrom_field"
// 1000G SNP reference panel directory: expects chr{N}.vcf.gz (+ .csi) per chromosome
params.snp_ref_panel_dir      = "/Users/thanhdng/Documents/databases/human/HC_1000G_hg38/rsid_vcfs"


// Other parameters
params.chromosomes           = 1..22
params.outdir                = "./imputation_output"
params.trace_dir             = "./imputation_traces"

// Imputation mode:
//   "2stage" (default) — Stage 1: SNP imputation per batch using 1000G panel,
//                        Stage 2: SV imputation on merged SNP-dense data
//   "1stage"           — SV imputation only (skip SNP pre-imputation)
//                        Useful when input data is already SNP-dense
params.imputation_mode       = "2stage"


log.info """\
================================================================
                        nf-imputation-sv
================================================================
                SV Imputation Pipeline with nextflow.
================================================================
"""

nextflow.enable.dsl=2

workflow {

    ref37  = file(params.genome_build_37)
    ref38  = file(params.genome_build_38)
    chain  = file(params.chain_file)

    bfile_ch = Channel
        .fromPath(params.input_manifest)
        .splitCsv(header: true)
        .filter { row -> row.prefix }
        .map { row ->
            def prefix = row.prefix
            def base   = prefix.tokenize('/')[-1]
            tuple(base, file("${prefix}.bed"), file("${prefix}.bim"), file("${prefix}.fam"))
        }

    // Step 1: Convert PLINK bfiles to normalised VCF (GRCh37)
    vcf_ch  = ArrayBfile_to_VCF(bfile_ch, ref37)

    // Step 2: Lift array VCFs from GRCh37 to GRCh38 using bcftools +liftover
    lifted_ch = LiftoverVcf(vcf_ch, chain, ref37, ref38)
        .map { batch, vcf, csi, rejected -> tuple(batch, vcf, csi) }

    // Step 3: Collect union of array positions (CHROM/POS) across all batches
    array_pos_ch = CollectArrayPositions( lifted_ch.map { batch, vcf, csi -> vcf }.collect() )

    ref_panel = file(params.impute_ref_panel)

    // Per-chromosome genetic map files staged as proper process inputs
    // Files: plink.chr{N}.GRCh38.map (numeric chrom IDs, no "chr" prefix)
    map_files_ch = Channel.fromList(params.chromosomes.toList())
        .map { chrom -> tuple(chrom, file("${params.plink_map_dir}/plink.chr${chrom}.GRCh38.map")) }

    if ( params.imputation_mode == "1stage" ) {

        // ── 1-stage mode: SV imputation directly on array data ────────────────
        // Step 4: Merge all batches into a single VCF per chromosome
        //   (no SNP pre-imputation; just combine batches for SV imputation)
        split_array_ch = SplitArrayVcf( lifted_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )
        // emits: (chrom, batch, vcf, csi)

        merged_array_input_ch = split_array_ch
            .map { chrom, batch, vcf, csi -> tuple(chrom, vcf, csi) }
            .groupTuple(by: 0)
        merged_array_ch = MergeArrayBatches(merged_array_input_ch)
        // emits: (chrom, merged_vcf, merged_csi)

        // Step 5: Filter SV ref panel to array positions + all SVs
        reduced_ref_ch = FilterRefPanel(array_pos_ch, ref_panel)

        // Step 6: Split SV ref panel by chromosome
        split_ref_ch = SplitRefPanel( reduced_ref_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )

        // Step 7: Impute SVs per chromosome on the merged array VCF
        BeagleSVImpute( merged_array_ch.combine(split_ref_ch, by: 0).combine(map_files_ch, by: 0) )

    } else {

        // ── 2-stage mode (default): SNP imputation → SV imputation ────────────
        // Step 4: Split lifted array VCFs by chromosome
        split_array_ch = SplitArrayVcf( lifted_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )
        // emits: (chrom, batch, vcf, csi)

        // Step 5: Build per-chromosome SNP reference panel
        snp_raw_ref_ch = Channel.fromList(params.chromosomes.toList())
            .map { chrom ->
                tuple(chrom,
                      file("${params.snp_ref_panel_dir}/chr${chrom}.vcf.gz"),
                      file("${params.snp_ref_panel_dir}/chr${chrom}.vcf.gz.csi"))
            }
        snp_ref_ch = BuildSNPRefPanel( array_pos_ch.combine(snp_raw_ref_ch) )
        // emits: (chrom, snp_ref_vcf, snp_ref_csi)

        // Step 6: Impute common SNPs per batch × chromosome
        snp_imputed_ch = BeagleSNPImpute( split_array_ch.combine(snp_ref_ch, by: 0).combine(map_files_ch, by: 0) )
        // emits: (chrom, batch, vcf, csi)

        // Step 7: Merge all batches per chromosome after SNP imputation
        snp_merge_input_ch = snp_imputed_ch
            .map { chrom, batch, vcf, csi -> tuple(chrom, vcf, csi) }
            .groupTuple(by: 0)
        merged_snp_ch = MergeSNPBatches(snp_merge_input_ch)
        // emits: (chrom, merged_vcf, merged_csi)

        // Step 8: Collect positions from all stage-1 merged VCFs
        sv_positions_ch = CollectSNPPositions( merged_snp_ch.map { chrom, vcf, csi -> vcf }.collect() )

        // Step 9: Filter SV ref panel to SNP positions + all SVs
        reduced_ref_ch = FilterRefPanel(sv_positions_ch, ref_panel)

        // Step 10: Split SV ref panel by chromosome
        split_ref_ch = SplitRefPanel( reduced_ref_ch.combine( Channel.fromList(params.chromosomes.toList()) ) )

        // Step 11: Impute SVs per chromosome on the merged SNP-imputed VCF
        BeagleSVImpute( merged_snp_ch.combine(split_ref_ch, by: 0).combine(map_files_ch, by: 0) )

    }
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

process MergeArrayBatches {

    tag "chr${chrom}"
    publishDir "${params.outdir}/array_merged", mode: 'copy'

    input:
    tuple val(chrom), path(vcfs), path(csis)

    output:
    tuple val(chrom), path("chr${chrom}.array_merged.vcf.gz"), path("chr${chrom}.array_merged.vcf.gz.csi")

    cpus 4
    memory '8 GB'

    script:
    def vcf_list = (vcfs instanceof List ? vcfs : [vcfs]).join(' ')
    def n        = (vcfs instanceof List ? vcfs : [vcfs]).size()
    """
    if [ ${n} -eq 1 ]; then
        cp ${vcf_list} chr${chrom}.array_merged.vcf.gz
    else
        bcftools merge \\
            --no-version \\
            --output-type z \\
            --output chr${chrom}.array_merged.vcf.gz \\
            --threads ${task.cpus} \\
            ${vcf_list}
    fi

    bcftools index -f chr${chrom}.array_merged.vcf.gz
    """
}

process CollectArrayPositions {

    publishDir "${params.trace_dir}/array_positions", mode: 'link'

    input:
    path vcfs  // all lifted array VCFs collected into one list

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

process BuildSNPRefPanel {

    tag "chr${chrom}"
    publishDir "${params.outdir}/snp_ref_by_chr", mode: 'link'

    input:
    tuple path(array_positions), val(chrom), path(snp_vcf), path(snp_csi)

    output:
    tuple val(chrom), path("chr${chrom}.snp_ref.vcf.gz"), path("chr${chrom}.snp_ref.vcf.gz.csi")

    cpus 2
    memory '8 GB'

    script:
    """
    # Extract array positions for this chromosome as a targets file (1-based CHROM\tPOS)
    awk -v c="${chrom}" '\$1 == c {print \$1"\t"\$2}' ${array_positions} > array_targets.txt

    # Pass 1: variants with MAF >= 5% (minor allele frequency across all samples)
    bcftools view \\
        --min-af 0.05:minor \\
        --output-type z \\
        --output maf5.vcf.gz \\
        ${snp_vcf}
    bcftools index -f maf5.vcf.gz

    # Pass 2: array positions regardless of MAF (to preserve typed variants)
    bcftools view \\
        --targets-file array_targets.txt \\
        --output-type z \\
        --output array_snps.vcf.gz \\
        ${snp_vcf}
    bcftools index -f array_snps.vcf.gz

    # Union: concatenate both passes, sort, and remove duplicates
    bcftools concat \\
        --allow-overlaps \\
        maf5.vcf.gz array_snps.vcf.gz \\
    | bcftools sort \\
    | bcftools norm \\
        --rm-dup all \\
        -Oz -o chr${chrom}.snp_ref.vcf.gz

    bcftools index -f chr${chrom}.snp_ref.vcf.gz

    rm -f maf5.vcf.gz maf5.vcf.gz.csi array_snps.vcf.gz array_snps.vcf.gz.csi array_targets.txt
    """
}

process CollectSNPPositions {

    publishDir "${params.trace_dir}/snp_positions", mode: 'link'

    input:
    path vcfs  // all stage-1 merged VCFs (one per chromosome) collected into one list

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
    tuple path(ref_vcf), path(ref_tbi), val(chrom)

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

process BeagleSNPImpute {

    tag "${batch}_chr${chrom}"
    publishDir "${params.outdir}/snp_imputed_vcfs", mode: 'copy'

    input:
    tuple val(chrom), val(batch), path(vcf), path(csi), path(snp_ref_vcf), path(snp_ref_csi), path(map_file)

    output:
    tuple val(chrom), val(batch),
          path("${batch}.chr${chrom}.snp_imputed.vcf.gz"),
          path("${batch}.chr${chrom}.snp_imputed.vcf.gz.csi")

    cpus 4
    memory '16 GB'

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

process BeagleSVImpute {

    tag "chr${chrom}"
    publishDir "${params.outdir}/sv_imputed_vcfs", mode: 'copy'

    input:
    tuple val(chrom), path(vcf), path(csi), path(sv_ref_vcf), path(sv_ref_csi), path(map_file)

    output:
    tuple val(chrom),
          path("chr${chrom}.sv_imputed.vcf.gz"),
          path("chr${chrom}.sv_imputed.vcf.gz.csi")

    cpus 4
    memory '16 GB'

    script:
    def mem_gb = (task.memory.toGiga() - 1) as int
    """
    java -Xmx${mem_gb}g -jar ${projectDir}/bin/beagle.27Feb25.75f.jar \\
        gt=${vcf} \\
        ref=${sv_ref_vcf} \\
        map=${map_file} \\
        chrom=${chrom} \\
        out=chr${chrom}.sv_imputed \\
        nthreads=${task.cpus} \\
        ap=true \\
        gp=true

    bcftools index -f chr${chrom}.sv_imputed.vcf.gz
    """
}

process LiftoverVcf {

    publishDir "${params.outdir}/lifted_vcfs", mode: 'copy'

    input:
    tuple val(batch), path(vcf), path(tbi)
    path chain
    path ref37        // source genome (GRCh37)
    path ref38        // destination genome (GRCh38)

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

process ArrayBfile_to_VCF {

    publishDir "${params.outdir}/array_vcfs", mode: 'copy'

    input:
    tuple val(batch), path(bed), path(bim), path(fam)
    path ref

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

    # Rename PLINK numeric codes to Ensembl notation (23→X, 24→Y, 26→MT; PAR already merged into X via --merge-x)
    printf '23\tX\n24\tY\n26\tMT\n' > chr_rename.txt

    # Compress VCF, rename chromosomes, then remove unzipped source
    bgzip -c ${batch}.vcf \
    | bcftools annotate \\
        --rename-chrs chr_rename.txt \\
        -Oz -o ${batch}.renamed.vcf.gz
    bcftools index -f ${batch}.renamed.vcf.gz

    # Remove intermediate files to save space
    rm -f ${batch}.vcf ${batch}.mergedx.bed ${batch}.mergedx.bim \
          ${batch}.mergedx.fam ${batch}.mergedx.nosex ${batch}.mergedx.log

    # Normalize: swap REF/ALT where possible, exclude remaining mismatches, remove duplicates
    bcftools norm \\
        -c sx \\
        --rm-dup all \\
        -f ${ref} \\
        -Oz \\
        -o ${batch}.norm.vcf.gz \\
        ${batch}.renamed.vcf.gz

    # Remove renamed intermediate
    rm -f ${batch}.renamed.vcf.gz ${batch}.renamed.vcf.gz.csi

    # Index
    bcftools index -f ${batch}.norm.vcf.gz

    # Quick sanity check: print any remaining REF mismatches (should be none)
    # || true prevents grep exit code 1 (no match = good) from failing the process
    bcftools norm \\
        -c e \\
        -f ${ref} \\
        ${batch}.norm.vcf.gz 2>&1 | grep -i "mismatch" || true
    """
}