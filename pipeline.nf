#!/usr/bin/env nextflow

// This is a comment.
// vim: syntax=groovy
// -*- mode: groovy;-*-

// You should configure everything in the '*.nfconfig' file, not here.

inputFastqs = Channel
    .fromFilePairs(params.fastq_files_glob,size: -1)
    { file -> file.getBaseName() - ~/_n0[12]/ - ~/.fastq/ }
    .ifEmpty{ error "couldn't find those fastqs!" }

file("./output").mkdirs()
file("./reports").mkdirs()

//  preprocessing
//
//	call_variants 1 # Call Variants Round 1
//	extract_snps 1 # Round 1. Extracts snps AND indels, separately
//	filter_snps 1 # Round 1
//	filter_indels 1 # Round 1
//
//	do_bqsr 1 # Do BQSR Round 1
//	do_bqsr 2 # Do BQSR Round 2
//
//	analyze_covariates # Only Done Once
//	apply_bqsr # Only Done Once
//
//	call_variants 2 # Call Variants Round 2
//	extract_snps 2 # Round 2. Extracts snps AND indels, separately
//	filter_snps 2 # Round 2
//	filter_indels 2 # Round 2
//
//	parse_metrics
//	do_snpeff

process align_to_reference {
    input:
        set val(pair_id), file(reads) from inputFastqs
    output:
        set val(pair_id), file("aligned_reads.sam") into alignedSam
    script:
    """
    module load ${params.modules.BWA}
    bwa mem -M -R "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.sequencing_platform}\\tPM:${params.sequencing_machine}\\tSM:${pair_id}" \
        ${params.reference_prefix} ${reads[0]} ${reads[1]} > \
        aligned_reads.sam
    """
}

process sort_and_bamize {
    input:
        set val(pair_id), file(aligned_reads) from alignedSam
    output:
        set val(pair_id), file("aligned_reads.bam") \
            into alignedSortedBam
    script:
    """
    module load ${params.modules.PICARD}
    java -jar ${params.modules.PICARD_JAR} SortSam \
        INPUT=${aligned_reads} OUTPUT=aligned_reads.bam \
        SORT_ORDER=coordinate
    """
}

( alignedSortedBam_for_metrics, alignedSortedBam_for_duplicates ) =
   alignedSortedBam.into(2)

process collect_alignment_metrics {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(aligned_reads) \
            from alignedSortedBam_for_metrics
    output:
        set val(pair_id), 
            file("${pair_id}_alignment_metrics.txt"), \
            file("${pair_id}_insert_metrics.txt"), \
            file("${pair_id}_insert_size_histogram.pdf"), \
            file("${pair_id}_depth_out.txt") \
            into alignmentMetrics_output
    script:
    """
    module load ${params.modules.PICARD}
    module load ${params.modules.R}
    module load ${params.modules.SAMTOOLS}
    java -jar ${params.modules.PICARD_JAR} \
        CollectAlignmentSummaryMetrics R=${params.reference_prefix} \
        I=${aligned_reads} O=${pair_id}_alignment_metrics.txt
    java -jar ${params.modules.PICARD_JAR} \
        CollectInsertSizeMetrics \
        INPUT=${aligned_reads} OUTPUT=${pair_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${pair_id}_insert_size_histogram.pdf 
    samtools depth -a ${aligned_reads} > ${pair_id}_depth_out.txt
    """
}

process mark_duplicates {
    publishDir "./output", mode: "copy", pattern: "*.txt"
    input:
        set val(pair_id), file(aligned_reads) \
            from alignedSortedBam_for_duplicates 
    output:
        set val(pair_id), file("${pair_id}_reads_dedup.bam"), \
            file("${pair_id}_reads_dedup.bai") \
            into alignedDedupedBam
    script:
    """
    module load ${params.modules.PICARD}
    java -jar ${params.modules.PICARD_JAR} MarkDuplicates \
        INPUT=${aligned_reads} OUTPUT=${pair_id}_reads_dedup.bam \
        METRICS_FILE=${pair_id}_metrics_dedup.txt
    java -jar ${params.modules.PICARD_JAR} BuildBamIndex \
        INPUT=${pair_id}_reads_dedup.bam
    """
}

process realign_around_indels {
    input:
        set val(pair_id), file(reads_dedup), file(reads_dedup_index) \
            from alignedDedupedBam
    output:
        set val(pair_id), val(1), \
            file("${pair_id}_realigned_reads.bam") \
            into realignedBam_first_output
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T RealignerTargetCreator \
        -R ${params.reference_prefix} \
        -I ${reads_dedup} \
        -o ${pair_id}_realignment_targets.list
    java -jar ${params.modules.GATK_JAR} -T IndelRealigner \
        -R ${params.reference_prefix} \
        -I ${reads_dedup} \
        -targetIntervals ${pair_id}_realignment_targets.list \
        -o ${pair_id}_realigned_reads.bam
    """
}

realignedBam_for_recalibration        = Channel.create()
realignedBam_for_reporting            = Channel.create()
realignedBam_for_variant_calling      = Channel.create()
realignedBam_for_base_recalibration   = Channel.create()
realignedBam_for_actual_recalibration = Channel.create()

realignedBam_first_output
    .tap(realignedBam_for_variant_calling)
    .tap(realignedBam_for_recalibration)


realignedBam_for_recalibration
    .tap(realignedBam_for_base_recalibration)
    .tap(realignedBam_for_actual_recalibration) 

process call_variants {
    input:
        set val(pair_id), val(round), file(input_bam) \
            from realignedBam_for_variant_calling
    output:
        set val(pair_id), val(round), file("${pair_id}_variants.vcf")\
            into calledVariants
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T HaplotypeCaller \
        -R ${params.reference_prefix} \
        -I ${input_bam} \
        -o ${pair_id}_variants.vcf
    """
}
calledVariants_for_extracting = Channel.create()
calledVariants_for_filtering = Channel.create()
calledVariants_final = Channel.create()
calledVariants
    .tap(calledVariants_for_extracting)
    .tap(calledVariants_for_filtering)
    .filter{it[1] == 2}
    .tap(calledVariants_final)
//    .filter{it[1] == 1}
//    .tap(calledVariants_for_base_recalibration)
process extract_snps_and_indels {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(variants) \
            from calledVariants_for_extracting
    output:
        set val(pair_id), val(round), \
            file("${pair_id}_snps_round${round}.vcf"), \
            file("${pair_id}_indel_round${round}.vcf") \
            into extractedVariants
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T SelectVariants \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -selectType SNP \
        -o ${pair_id}_snps_round${round}.vcf
    java -jar ${params.modules.GATK_JAR} -T SelectVariants \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -selectType INDEL \
        -o ${pair_id}_indel_round${round}.vcf
    """
}
process filter_snps_and_indels {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(variants) \
            from calledVariants_for_filtering
    output:
        set val(pair_id), val(round), 
            file("${pair_id}_filtered_snps_round${round}.vcf"), \
            file("${pair_id}_filtered_indels_round${round}.vcf") \
            into filteredVariants
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T VariantFiltration \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -filterName \"QD_filter\" \
        -filter \"QD < 2.0\" \
        -filterName \"FS_filter\" \
        -filter \"FS > 60.0\" \
        -filterName \"MQ_filter\" \
        -filter \"MQ < 40.0\" \
        -filterName \"SOR_filter\" \
        -filter \"SOR > 4.0\" \
        -o ${pair_id}_filtered_snps_round${round}.vcf
    java -jar ${params.modules.GATK_JAR} -T VariantFiltration \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -filterName \"QD_filter\" \
        -filter \"QD < 2.0\" \
        -filterName \"FS_filter\" \
        -filter \"FS > 200.0\" \
        -filterName \"SOR_filter\" \
        -filter \"SOR > 10.0\" \
        -o ${pair_id}_filtered_indels_round${round}.vcf
    """
}
filteredVariants_for_recalibration = Channel.create()
filteredVariants_for_final_report  = Channel.create()
filteredVariants
    .choice(filteredVariants_for_recalibration,
        filteredVariants_for_final_report)
    { a -> a[1] > 1 ? 0 : 1 }

process base_recalibrator {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(input_bam) \
            from realignedBam_for_base_recalibration
        set val(pair_id), val(round), \
            file(filtered_snps), file(filtered_indels) \
            from filteredVariants_for_recalibration
    output:
        set val(pair_id), file("first_recal_data.table"), \
            file("second_recal_data.table") \
            into bqsrOutputs
//	#todo: knownSites input shouldnt be full raw_variants.vcf file but only the TOP variants!
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T BaseRecalibrator \
        -R ${params.reference_prefix} \
        -I ${input_bam} \
        -knownSites ${filtered_snps} -knownSites ${filtered_indels} \
        -o first_recal_data.table
    java -jar ${params.modules.GATK_JAR} -T BaseRecalibrator \
        -R ${params.reference_prefix} \
        -I ${input_bam} \
        -knownSites ${filtered_snps} -knownSites ${filtered_indels} \
        -BQSR first_recal_data.table
        -o second_recal_data.table
    """
}
bqsrOutputs_for_covariate = Channel.create()
bqsrOutputs_for_apply = Channel.create()
bqsrOutputs
    .tap(bqsrOutputs_for_covariate)
    .tap(bqsrOutputs_for_apply)

process covariates_analyzer {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(first_recal_table), \
            file(second_recal_table) \
            from bqsrOutputs_for_covariate
    output:
        set val(pair_id), file("${pair_id}_recalibration_plots.pdf") \
            into analyzedCovariates_output 
    script:
    """
    module load ${params.modules.GATK}
    module load ${params.modules.R}
    java -jar ${params.modules.GATK_JAR} -T AnalyzeCovariates \
        -R ${params.reference_prefix} \
        -before ${first_recal_table} \
        -after ${second_recal_table} \
        -plots ${pair_id}_recalibration_plots.pdf
    """
}

process apply_bqsr {
    input:
        set val(pair_id), \
            file(first_recal_table), file(second_recal_table) \
            from bqsrOutputs_for_apply
        set val(pair_id), val(round), file(input_bam) \
            from realignedBam_for_actual_recalibration
    output:
        set val(pair_id), val(2), \
            file("${pair_id}_recalibrated.bam") \
            into realignedBam_second_pass
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T PrintReads \
        -R ${params.reference_prefix} \
        -I ${input_bam} \
        -BQSR ${first_recal_table} \
        -o ${pair_id}_recalibrated.bam
    """
}

realignedBam_second_pass

process call_variants2 {
    input:
        set val(pair_id), val(round), file(input_bam) \
            from realignedBam_second_pass
    output:
        set val(pair_id), val(round), file("${pair_id}_variants.vcf")\
            into calledVariants2
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T HaplotypeCaller \
        -R ${params.reference_prefix} \
        -I ${input_bam} \
        -o ${pair_id}_variants.vcf
    """
}
calledVariants2_for_extracting = Channel.create()
calledVariants2_for_filtering = Channel.create()
calledVariants2_final = Channel.create()
calledVariants2
    .tap(calledVariants2_for_extracting)
    .tap(calledVariants2_for_filtering)
    .filter{it[1] == 2}
    .tap(calledVariants2_final)
//    .filter{it[1] == 1}
//    .tap(calledVariants_for_base_recalibration)
process extract_snps_and_indels2 {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(variants) \
            from calledVariants2_for_extracting
    output:
        set val(pair_id), val(round), \
            file("${pair_id}_snps_round${round}.vcf"), \
            file("${pair_id}_indel_round${round}.vcf") \
            into extractedVariants2
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T SelectVariants \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -selectType SNP \
        -o ${pair_id}_snps_round${round}.vcf
    java -jar ${params.modules.GATK_JAR} -T SelectVariants \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -selectType INDEL \
        -o ${pair_id}_indel_round${round}.vcf
    """
}
process filter_snps_and_indels2 {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(variants) \
            from calledVariants2_for_filtering
    output:
        set val(pair_id), val(round), 
            file("${pair_id}_filtered_snps_round${round}.vcf"), \
            file("${pair_id}_filtered_indels_round${round}.vcf") \
            into filteredVariants2
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T VariantFiltration \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -filterName \"QD_filter\" \
        -filter \"QD < 2.0\" \
        -filterName \"FS_filter\" \
        -filter \"FS > 60.0\" \
        -filterName \"MQ_filter\" \
        -filter \"MQ < 40.0\" \
        -filterName \"SOR_filter\" \
        -filter \"SOR > 4.0\" \
        -o ${pair_id}_filtered_snps_round${round}.vcf
    java -jar ${params.modules.GATK_JAR} -T VariantFiltration \
        -R ${params.reference_prefix} \
        -V ${variants} \
        -filterName \"QD_filter\" \
        -filter \"QD < 2.0\" \
        -filterName \"FS_filter\" \
        -filter \"FS > 200.0\" \
        -filterName \"SOR_filter\" \
        -filter \"SOR > 10.0\" \
        -o ${pair_id}_filtered_indels_round${round}.vcf
    """
}
filteredVariants2_for_recalibration = Channel.create()
filteredVariants2_for_final_report  = Channel.create()
filteredVariants2
    .choice(filteredVariants2_for_recalibration,
        filteredVariants2_for_final_report)
    { a -> a[1] > 1 ? 0 : 1 }

process final_metrics {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(recalibrated_bam) \
            from realignedBam_for_reporting
    output:
        set val(pair_id), file("${pair_id}_genomecov.bedgraph") \
            into final_metrics_output
    script:
    """
    module load ${params.modules.BEDTOOLS}
    bedtools genomecov -bga -ibam ${recalibrated_bam} > \
        ${pair_id}_genomecov.bedgraph 
    #sh ./parse_metrics.sh ${pair_id}"
    """
}

process snpeff {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(final_variants) \
            from calledVariants_final
    output:
        set val(pair_id), val(round),
            file("${pair_id}_filtered_snps_final.ann.vcf") \
            into calledVariants_output
    script:
    """
    module load ${params.modules.SNPEFF}
    java -jar ${params.modules.SNPEFF_JAR} \
        -v ${params.snpeff_database} ${final_variants} > \
        ${pair_id}_filtered_snps_final.ann.vcf"
    """
}

