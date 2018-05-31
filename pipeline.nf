#!/usr/bin/env nextflow

// This is a comment.
// vim: syntax=groovy
// -*- mode: groovy;-*-

// You should configure everything in the '*.nfconfig' file, not here.

input_fastqs = Channel
    .fromFilePairs(params.fastq_files_glob,size: -1)
    { file -> file.getBaseName() - ~/_n0[12]/ - ~/.fastq/ }
    .ifEmpty{ error "couldn't find those fastqs!" }

file("./output").mkdirs()
file("./reports").mkdirs()


//  preprocessing
//	call_variants 1 # Call Variants Round 1
//	extract_snps 1 # Round 1. Extracts snps AND indels, separately
//	filter_snps 1 # Round 1
//	filter_indels 1 # Round 1
//	do_bqsr 1 # Do BQSR Round 1
//	do_bqsr 2 # Do BQSR Round 2
//	analyze_covariates # Only Done Once
//	apply_bqsr # Only Done Once
//	call_variants 2 # Call Variants Round 2
//	extract_snps 2 # Round 2. Extracts snps AND indels, separately
//	filter_snps 2 # Round 2
//	filter_indels 2 # Round 2
//	parse_metrics
//	do_snpeff
//done

process align_to_reference {
    input:
        set val(pair_id), file(reads) from input_fastqs
    output:
        set val(pair_id), file("aligned_reads.sam") into aligned_reads
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
        set val(pair_id), file(aligned_reads) from aligned_reads
    output:
        set val(pair_id), file("aligned_reads.bam") \
            into aligned_reads_sorted
    script:
    """
    module load ${params.modules.PICARD}
    java -jar ${params.modules.PICARD_JAR} SortSam \
        INPUT=${aligned_reads} OUTPUT=aligned_reads.bam \
        SORT_ORDER=coordinate
    """
}

( aligned_reads_sorted_metrics, aligned_reads_sorted_duplicates ) =
   aligned_reads_sorted.into(2)

process collect_alignment_metrics {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(aligned_reads) \
            from aligned_reads_sorted_metrics
    output:
        set val(pair_id), 
            file("${pair_id}_alignment_metrics.txt"), \
            file("${pair_id}_insert_metrics.txt"), \
            file("${pair_id}_insert_size_histogram.pdf"), \
            file("${pair_id}_depth_out.txt") \
            into alignment_metrics
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
            from aligned_reads_sorted_duplicates 
    output:
        set val(pair_id), file("${pair_id}_reads_dedup.bam"), \
            file("${pair_id}_reads_dedup.bai") \
            into aligned_reads_deduplicated
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
            from aligned_reads_deduplicated
    output:
        set val(pair_id), val(1), \
            file("${pair_id}_realigned_reads.bam") \
            into realigned_around_indels
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

//input_to_call_variants = Channel.create()
//downstream = Channel.create()
//realigned_around_indels
//    .choice(input_to_call_variants,downstream)
//    { a -> a[1] == 1 ? 0 : 1 }

( input_to_call_variants, input_to_base_recalibration, 
    input_to_actual_recalibration ) =
    realigned_around_indels.into(3)

process call_variants {
    input:
        set val(pair_id), val(round), file(input_bam) \
            from input_to_call_variants
    output:
        set val(pair_id), val(round), file("${pair_id}_variants.vcf")\
            into variants_called
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T HaplotypeCaller \
        -R ${params.reference_prefix} \
        -I ${input_bam} \
        -o ${pair_id}_variants.vcf
    """
}

( variants_called_for_extracting, variants_called_for_filtering) = 
   variants_called.into(2)

process extract_snps_and_indels {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(variants) \
            from variants_called_for_extracting
    output:
        set val(pair_id), \
            file("${pair_id}_snps_round${round}.vcf"), \
            file("${pair_id}_indel_round${round}.vcf") \
            into variants_extracted
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
            from variants_called_for_filtering
    output:
        set val(pair_id), 
            file("${pair_id}_filtered_snps_round${round}.vcf"), \
            file("${pair_id}_filtered_indels_round${round}.vcf") \
            into variants_filtered
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

process base_recalibrator {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(input_bam) \
            from input_to_base_recalibration
        set val(pair_id), file(filtered_snps), file(filtered_indels) \
            from variants_filtered
    output:
        set val(pair_id), file("first_recal_data.table"),
            file("second_recal_data.table")
            into bqsr_recalibrated
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

(bqsr_recalibrated_covariate, bqsr_recalibrated_apply ) =
    bqsr_recalibrated.into(2)

process covariates_analyzer {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(first_recal_table),
            file(second_recal_table)
            from bqsr_recalibrated_covariate
    output:
        set val(pair_id), file("${pair_id}_recalibration_plots.pdf"),
            into covariates_analyzed
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
        set val(pair_id), file(first_recal_table),
            file(second_recal_table)
            from bqsr_recalibrated_apply
        set val(pair_id), val(round), file(input_bam) \
            from input_to_actual_recalibration
    output:
        set val(pair_id), file("${pair_id}_recalibrated.bam") \
            into bqsr_recalibrated_reads
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

//parse_metrics(){
//com="cd $FWD && \
//module load $BEDTOOLS && \
//bedtools genomecov -bga -ibam ${ID}_recal_reads.bam > ${ID}_genomecov.bedgraph && \
//sh /scratch/work/cgsb/scripts/variant_calling/parse_metrics.sh $ID"
//response=\
//$(sbatch -J $ID.parseMetrics -o $ID.parseMetrics.out -e $ID.parseMetrics.err --dependency=afterok:$filterSnps_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
//stringarray=($response)
//parseMetrics=${stringarray[-1]}
//echo $parseMetrics >> $ID.log
//echo "PARSEMETRICS: " $parseMetrics
//echo "ParseMetrics Submitted"
//}

//do_snpeff(){
//com="cd $FWD && \
//module load $SNPEFF && \
//java -jar $SNPEFF_JAR \
//-v $SNPEFF_DB \
//${ID}_filtered_snps_final.vcf > ${ID}_filtered_snps_final.ann.vcf"
//response=\
//$(sbatch -J $ID.snpEff -o $ID.snpEff.out -e $ID.snpEff.err --dependency=afterok:$filterSnps_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
//stringarray=($response)
//snpEff=${stringarray[-1]}
//echo $snpEff >> $ID.log
//echo "SNPEFF: " $snpEff
//echo "SnpEFF submitted"
//}
//
//## Build Report Header and Create File
//REPORT_HEADER="ID,# reads,aligned reads,% aligned,aligned bases,read length,% paired,mean insert size,# SNPs 1,# SNPs 1 filtered,# SNPs 2, # SNPs filtered 2,average coverage"
//echo $REPORT_HEADER > report.csv
//
//
