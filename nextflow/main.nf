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
    convert_snps --input ${snps} --output variants.vcf --reference ${reference}
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
    container params.containers.bcftools
    
    input:
    path vcf
    path reference
    
    output:
    path "normalized.vcf.gz", emit: normalized_vcf
    path "normalized.vcf.gz.tbi", emit: normalized_vcf_index
    
    script:
    """
    # Get reference name from FASTA header
    REF_NAME=\$(head -n1 ${reference} | sed 's/^>//' | cut -f1 -d' ')
    REF_LENGTH=\$(grep -v "^>" ${reference} | tr -d '\n' | wc -c)
    
    # Create a VCF header with proper contig definition
    cat > header.txt << EOH
##fileformat=VCFv4.2
##reference=MG1655.fna
##contig=<ID=\${REF_NAME},length=\${REF_LENGTH}>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOH

    # Fix contig name in VCF and combine with new header
    grep -v "^#" ${vcf} | sed "s/MG1655.fna/\${REF_NAME}/" > variants.txt
    
    # Create header with proper tab characters
    { 
        cat header.txt
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        cat variants.txt
    } > fixed.vcf
    
    # Sort and normalize
    bcftools sort fixed.vcf > sorted.vcf
    bcftools norm -f ${reference} sorted.vcf > normalized.vcf
    bgzip -f normalized.vcf
    tabix -p vcf normalized.vcf.gz
    """
}

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
