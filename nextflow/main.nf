#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.reference = null
params.query = null
params.outdir = "results"
params.aligner = "mummer"  // options: "mummer" or "minimap2"
params.min_sv_size = 50
params.max_sv_size = 100000

// Parameter validation
if (!params.reference || !params.query) {
    error "Please provide reference and query genome files using --reference and --query"
}
if (!params.aligner in ['mummer', 'minimap2']) {
    error "Aligner must be either 'mummer' or 'minimap2'"
}

def helpMessage() {
    log.info"""
    Variant Calling Pipeline
    =======================
    A Nextflow pipeline for variant calling using either MUMmer or Minimap2.

    Required Arguments:
    --reference      Path to reference genome (FASTA)
    --query         Path to query genome (FASTA)

    Optional Arguments:
    --aligner       Alignment tool to use (mummer/minimap2) [default: mummer]
    --outdir        Output directory [default: results]
    --detect_sv     Enable structural variant detection [default: false]
    
    SV Detection Options (only used when --detect_sv is enabled):
    --min_sv_size   Minimum size for structural variants [default: 50]
    --max_sv_size   Maximum size for structural variants [default: 100000]

    Example Usage:
    # Basic SNP calling with MUMmer
    nextflow run main.nf --reference ref.fa --query query.fa

    # Include structural variant detection
    nextflow run main.nf --reference ref.fa --query query.fa --detect_sv

    # Use minimap2 instead of MUMmer
    nextflow run main.nf --reference ref.fa --query query.fa --aligner minimap2
    """
}

// Parameter validation
def validateParameters() {
    // Show help if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Check required parameters
    if (!params.reference || !params.query) {
        log.error "Please provide reference and query genome files using --reference and --query"
        exit 1
    }

    // Validate aligner choice
    if (!params.aligner in ['mummer', 'minimap2']) {
        log.error "Aligner must be either 'mummer' or 'minimap2'"
        exit 1
    }

    // Validate SV parameters if SV detection is enabled
    if (params.detect_sv) {
        if (params.min_sv_size <= 0) {
            log.error "Minimum SV size must be greater than 0"
            exit 1
        }
        if (params.max_sv_size <= params.min_sv_size) {
            log.error "Maximum SV size must be greater than minimum SV size"
            exit 1
        }
        log.info "SV detection enabled (size range: ${params.min_sv_size} - ${params.max_sv_size}bp)"
    }
}

// MUMmer Processes
process NUCMER_ALIGN {
    publishDir "${params.outdir}/nucmer", mode: 'copy'
    container params.containers.mummer
    
    input:
    path reference
    path query
    
    output:
    path "output.delta", emit: delta
    
    script:
    """
    nucmer --maxmatch -p output ${reference} ${query}
    """
}

process FILTER_DELTA {
    publishDir "${params.outdir}/filtered", mode: 'copy'
    container params.containers.mummer
    
    input:
    path delta
    
    output:
    path "filtered.delta", emit: filtered_delta
    
    script:
    """
    delta-filter -r -q ${delta} > filtered.delta
    """
}

process SHOW_SNPS {
    publishDir "${params.outdir}/snps", mode: 'copy'
    container params.containers.mummer
    
    input:
    path filtered_delta
    path reference
    
    output:
    path "variants.snps", emit: snps
    
    script:
    """
    show-snps -Clr -T ${filtered_delta} > variants.snps
    """
}

process SHOW_COORDS {
    publishDir "${params.outdir}/coords", mode: 'copy'
    container params.containers.mummer
    
    input:
    path filtered_delta
    
    output:
    path "coords.txt", emit: coords
    
    script:
    """
    show-coords -THrd ${filtered_delta} > coords.txt
    """
}

process MINIMAP_ALIGN {
    publishDir "${params.outdir}/minimap2", mode: 'copy'
    container params.containers.minimap2
    
    input:
    path reference
    path query
    
    output:
    path "alignment.paf", emit: paf
    
    script:
    """
    minimap2 -cx asm5 --cs -t ${task.cpus} ${reference} ${query} > alignment.paf
    """
}

process PAF_TO_VCF {
    publishDir "${params.outdir}/minimap2", mode: 'copy'
    container params.containers.variant_caller
    
    input:
    path paf
    path reference
    
    output:
    path "minimap2_variants.vcf", emit: vcf
    
    script:
    """
    python -m variant_caller.paf_to_vcf \
        --input ${paf} \
        --output minimap2_variants.vcf \
        --reference ${reference}
    """
}

process DETECT_STRUCTURAL_VARIANTS {
    publishDir "${params.outdir}/sv", mode: 'copy'
    container params.containers.variant_caller
    
    input:
    path coords
    
    output:
    path "structural_variants.txt", emit: svs
    
    script:
    """
    detect_sv --input ${coords} --output structural_variants.txt \
        --min-size ${params.min_sv_size} --max-size ${params.max_sv_size}
    """
}

process CONVERT_TO_VCF {
    publishDir "${params.outdir}/vcf", mode: 'copy'
    container params.containers.variant_caller
    
    input:
    path snps
    path reference
    
    output:
    path "variants.vcf", emit: vcf
    
    script:
    """
    #!/usr/bin/env python3
    from variant_caller.vcf_converter import convert_snps_to_vcf
    
    convert_snps_to_vcf("${snps}", "${reference}", "variants.vcf")
    """
}

process MERGE_VARIANTS_TO_VCF {
    publishDir "${params.outdir}/final", mode: 'copy'
    container params.containers.variant_caller
    
    input:
    path snps_vcf
    path svs
    path reference
    
    output:
    path "merged_variants.vcf", emit: merged_vcf
    
    script:
    """
    merge_vcf --snps ${snps_vcf} --svs ${svs} --output merged_variants.vcf --reference ${reference}
    """
}

process SORT_AND_NORMALIZE_VCF {
    publishDir "${params.outdir}/final", mode: 'copy'
    container params.containers.variant_caller  // Changed from bcftools to variant_caller
    
    input:
    path vcf
    path reference
    
    output:
    path "normalized.vcf.gz", emit: normalized_vcf
    path "normalized.vcf.gz.tbi", emit: normalized_vcf_index
    
    script:
    """
    #!/usr/bin/env python3
    
    from variant_caller.vcf_utils import VCFHandler
    import logging
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    try:
        # Initialize VCF handler
        handler = VCFHandler("${reference}")
        
        # Sort and normalize VCF
        logger.info("Sorting and normalizing VCF file...")
        handler.sort_and_normalize("${vcf}", "normalized.vcf")
        
        # Compress output
        logger.info("Compressing output...")
        import subprocess
        subprocess.run(["bgzip", "-f", "normalized.vcf"], check=True)
        subprocess.run(["tabix", "-p", "vcf", "normalized.vcf.gz"], check=True)
        
    except Exception as e:
        logger.error(f"Error processing VCF: {e}")
        raise
    """
}

// Workflow definitions
workflow minimap_flow {
    take:
        reference
        query
    
    main:
        MINIMAP_ALIGN(reference, query)
        PAF_TO_VCF(MINIMAP_ALIGN.out.paf, reference)
        SORT_AND_NORMALIZE_VCF(PAF_TO_VCF.out.vcf, reference)
        
    emit:
        vcf = SORT_AND_NORMALIZE_VCF.out.normalized_vcf
        index = SORT_AND_NORMALIZE_VCF.out.normalized_vcf_index
}

workflow mummer_flow {
    take:
        reference
        query
    
    main:
        NUCMER_ALIGN(reference, query)
        FILTER_DELTA(NUCMER_ALIGN.out.delta)
        SHOW_SNPS(FILTER_DELTA.out.filtered_delta, reference)
        CONVERT_TO_VCF(SHOW_SNPS.out.snps, reference)
        SHOW_COORDS(FILTER_DELTA.out.filtered_delta)
        DETECT_STRUCTURAL_VARIANTS(SHOW_COORDS.out.coords)
        MERGE_VARIANTS_TO_VCF(
            CONVERT_TO_VCF.out.vcf,
            DETECT_STRUCTURAL_VARIANTS.out.svs,
            reference
        )
        SORT_AND_NORMALIZE_VCF(MERGE_VARIANTS_TO_VCF.out.merged_vcf, reference)
    
    emit:
        vcf = SORT_AND_NORMALIZE_VCF.out.normalized_vcf
        index = SORT_AND_NORMALIZE_VCF.out.normalized_vcf_index
}

workflow {
    reference_ch = channel.fromPath(params.reference)
    query_ch = channel.fromPath(params.query)
    
    if (params.aligner == 'mummer') {
        mummer_flow(reference_ch, query_ch)
    } else {
        minimap_flow(reference_ch, query_ch)
    }
}
