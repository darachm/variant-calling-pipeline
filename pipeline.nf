#!/usr/bin/env nextflow

// This is a comment. I like to comment before the line it refers to.

// You should configure everything in the '*.nfconfig' file, not here.
// TODO put all config in this file?

// This is first a channel we're going to use later to count the 
// inputs. More details later.
input_count_channel = Channel.create()

// This block uses the 'params.fastq_files_glob' to find fastq files.
// It then uses the size parameter to not auto-group, and instead uses
// the mapping through '.getBaseName()' and subtracts two regexs to
// trim that name down to an ID string. Then it uses 'tap' to copy
// it into the 'input_count_channel'. '.ifEmpty', it reports an error.
input_fastqs = Channel
    .fromFilePairs( params.fastq_files_glob , size: -1)
    { file -> file.getBaseName() - ~/_n0[12]/ - ~/.fastq/ }
    .tap(input_count_channel)
    .ifEmpty{ error "couldn't find those fastqs!" }

// Now, we make an variable, then 'subscribe' to the input count
// channel. This consumes the channel, so that's why we're using a
// 'tap'ed copy. Everytime something comes through, we increment it
// by one to count it.
input_count = 0
input_count_channel.subscribe({ input_count += 1 })

// These are two directories for outputs and reports.
file("./output").mkdirs()
file("./reports").mkdirs()

// The first process takes the raw input fastqs and aligns them to
// the reference to make a sam file.
// Note that the ID used to group the paired-end reads is kept in the
// first position. Because of this, most of the channels are 'set's
// of 'val'ues and 'file's. 
process align_to_reference {
    input:
        set val(pair_id), file(reads) from input_fastqs
    output:
        set val(pair_id), file("aligned_reads.sam") into aligned_sam
    script:
    """
    module load ${params.modules.BWA}
    bwa mem -M -R "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.sequencing_platform}\\tPM:${params.sequencing_machine}\\tSM:${pair_id}" \
        ${params.reference_prefix} ${reads[0]} ${reads[1]} > \
        aligned_reads.sam
    """
}

// The sam file is sorted into a bam file.
process sort_and_bamize {
    input:
        set val(pair_id), file(aligned_reads) from aligned_sam
    output:
        set val(pair_id), file("aligned_reads.bam") \
            into aligned_sorted_bam
    script:
    """
    module load ${params.modules.PICARD}
    java -jar ${params.modules.PICARD_JAR} SortSam \
        INPUT=${aligned_reads} OUTPUT=aligned_reads.bam \
        SORT_ORDER=coordinate
    """
}

// The channel of sorted bam files is split into two channels, one for
// calculating metrics, and one for deduplicating. Note that both are
// channels of 'set's of input ID, round ID, and a bam filepath.
( aligned_sorted_bam_for_metrics, aligned_sorted_bam_for_duplicates ) =
   aligned_sorted_bam.into(2)

// On the metrics copy, calculate the alignment metrics. 
// Note this is copied ot the output directory
process collect_alignment_metrics {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(aligned_reads) \
            from aligned_sorted_bam_for_metrics
    output:
        set val(pair_id), 
            file("${pair_id}_alignment_metrics.txt"), \
            file("${pair_id}_insert_metrics.txt"), \
            file("${pair_id}_insert_size_histogram.pdf"), \
            file("${pair_id}_depth_out.txt") \
            into alignment_metrics_output
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

// On the deduplicate copy, do the deduplication and copy the output
// txt files to the output directory.
process mark_duplicates {
    publishDir "./output", mode: "copy", pattern: "*.txt"
    input:
        set val(pair_id), file(aligned_reads) \
            from aligned_sorted_bam_for_duplicates 
    output:
        set val(pair_id), file("${pair_id}_reads_dedup.bam"), \
            file("${pair_id}_reads_dedup.bai") \
            into aligned_deduped_bam
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

// Take that deduplicated file and realign it around the indels.
process realign_around_indels {
    input:
        set val(pair_id), file(reads_dedup), file(reads_dedup_index) \
            from aligned_deduped_bam 
    output:
        set val(pair_id), val(1), \
            file("${pair_id}_realigned_reads.bam") \
            into bams_to_work_on_first_input
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

// Now it gets a bit tricky. We want to take this first rounds of
// inputs and process them, but also leave it open to 'mix'ing in
// some downstream files in a bit of a loop.

// First, we make all the channels. Note the second one is a throwaway
// '*_tmp' channel, because nextflow is careful about single 
// inputs/outputs.
bams_for_variant_calling = Channel.create()
bams_for_variant_calling_tmp = Channel.create()
bams_for_base_recalibration = Channel.create()
bams_for_actual_recalibration = Channel.create()
bam_recalibrated_qualitites  = Channel.create()

// Here we actually 'mix' in the bams with recalibrated quality 
// scores.  Note that everything is a 'set' with the second position
// (as 1, they are 0-indexed) value as the rounds (1-indexed for the
// end-users). This always 'tap's the channel as a copy into a bam
// for variant calling and a dummy copy for something later.
// Then, it runs a closure in the curly braces. If we are on round 1,
// so the current set called 'it' has the second value 'it[1]' equal
// to 1, then this is filtered into the bams for recalibrations.
bams_to_work_on_first_input.mix(bam_recalibrated_qualitites)
    .tap(bams_for_variant_calling)
    .tap(bams_for_variant_calling_tmp)
    .filter({ it[1] == 1 })
    .tap(bams_for_base_recalibration)
    .tap(bams_for_actual_recalibration)

// Here we used that dummy variable to make a bam channel for final
// bedgraph generation, if we are on round 2. Everything a channel is
// used, it is consumed, so this is why we need a dummy copy before
// filtering it.
bams_for_reporting = Channel.create()
bams_for_variant_calling_tmp
    .filter({ it[1] == 2 })
    .tap(bams_for_reporting)

// On everything, we call variants.
process call_variants {
    input:
        set val(pair_id), val(round), file(input_bam) \
            from bams_for_variant_calling
    output:
        set val(pair_id), val(round), file("${pair_id}_variants.vcf")\
            into called_variants
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T HaplotypeCaller \
        -R ${params.reference_prefix} \
        -I ${input_bam} \
        -o ${pair_id}_variants.vcf
    """
}

// We then do a similar thing as before. We make a few channels,
// then 'tap' into them, then 'filter', then 'tap' if it passes that
// 'filter'. This way, the variants are called and filtered on both
// rounds, but only the final round has the variant files for snpeff.
called_variants_for_extracting = Channel.create()
called_variants_for_filtering = Channel.create()
called_variants_for_reporting = Channel.create()
called_variants
    .tap(called_variants_for_extracting)
    .tap(called_variants_for_filtering)
    .filter({ it[1] == 2 })
    .tap(called_variants_for_reporting)

// Extract snp and indel variants, publish the files by copying to
// output directory.
process extract_snps_and_indels {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(variants) \
            from called_variants_for_extracting
    output:
        set val(pair_id), val(round), \
            file("${pair_id}_snps_round${round}.vcf"), \
            file("${pair_id}_indel_round${round}.vcf") \
            into extracted_variants
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

// And this filters them. It also outputs a channel of the filtered
// variants for use by the recalibration steps.
process filter_snps_and_indels {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(variants) \
            from called_variants_for_filtering
    output:
        set val(pair_id), val(round), 
            file("${pair_id}_filtered_snps_round${round}.vcf"), \
            file("${pair_id}_filtered_indels_round${round}.vcf") \
            into filtered_variants
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

// Now, we make a copy of these filtered variants that were done on
// the first round. 
filtered_variants_for_recalibration = filtered_variants
    .filter({it[1] == 1})
// Then we join the bams allocated for the recalibration with these
// filtered first round of variants. By default it joins on the value
// in the first position of the set, but here I've just done that
// same thing explicity by setting the 'by' parameter. 
// NOTE the colon, not =.
bam_and_variants_for_recalibration = bams_for_base_recalibration
    .join(filtered_variants_for_recalibration, by: 0)

// Then we can apply the base recalibrator.
process base_recalibrator {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(input_bam), \
            val(round2), file(filtered_snps), file(filtered_indels) \
            from bam_and_variants_for_recalibration
    output:
        set val(pair_id), file("first_recal_data.table"), \
            file("second_recal_data.table") \
            into bqsr_outputs 
// todo from the bash script: 
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
        -BQSR first_recal_data.table \
        -o second_recal_data.table
    """
}

// Then we need to split ths output because it's consumed each time.
( bqsr_outputs_for_covariate, bqsr_outputs_for_apply 
    ) = bqsr_outputs.into(2)

// Here's the covariates analysis.
process covariates_analyzer {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(first_recal_table), \
            file(second_recal_table) \
            from bqsr_outputs_for_covariate
    output:
        set val(pair_id), file("${pair_id}_recalibration_plots.pdf") \
            into analyzed_covariates_output 
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

// Then to apply the actual recalibration we join the bams allocated
// to the bqsr outputs, and that's joined on that ID variable in the
// first position.
bam_and_bqsr_for_actual_recalibration = bams_for_actual_recalibration
    .join(bqsr_outputs_for_apply)

// Then we apply the recalibration.
process apply_bqsr {
    input:
        set val(pair_id), val(round), file(input_bam), \
            file(first_recal_table), file(second_recal_table) \
            from bam_and_bqsr_for_actual_recalibration
    output:
        set val(pair_id), val(2), \
            file("${pair_id}_recalibrated.bam") \
            into bam_recalibrated_qualitites
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

// This still needs to be debugged, and the parse_metrics script 
// re-writen.
//process final_metrics {
//    publishDir "./output", mode: "copy"
//    input:
//        set val(pair_id), val(round), file(recalibrated_bam) \
//            from bams_for_reporting 
//    output:
//        set val(pair_id), file("${pair_id}_round${round}_genomecov.bedgraph") \
//            into output_bedgraph
//    script:
//    """
//    module load ${params.modules.BEDTOOLS}
//    bedtools genomecov -bga -ibam ${recalibrated_bam} > \
//        ${pair_id}_round${round}_ genomecov.bedgraph 
//    #sh ./parse_metrics.sh ${pair_id}"
//    """
//}

// This runs the SNPEFF step, waiting for variants from the second
// round before running.
process snpeff {
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), val(round), file(final_variants) \
            from called_variants_for_reporting
    output:
        set val(pair_id), val(round),
            file("${pair_id}_filtered_snps_final.ann.vcf") \
            into called_variants_output
    script:
    """
    module load ${params.modules.SNPEFF}
    java -jar ${params.modules.SNPEFF_JAR} \
        -v ${params.snpeff_database} ${final_variants} > \
        ${pair_id}_filtered_snps_final.ann.vcf
    """
}

// This part is essential for the looping. When we used 'mix' above,
// nextflow left the channel open. Here, we count the outputs by
// using subscribe. In that same closure, then we test if the output
// count is equal to the input count. If so, then we send the close()
// signal to the channel that's left open. Because we also 'tap'ed 
// that to other channels, we also need to close those. Then it
// should hang a bit, then close when it's decided that the flow is
// over.
output_count = 0
called_variants_output.subscribe({ 
        output_count += 1 ; 
        if (output_count == input_count) {
            bams_to_work_on_first_input.close()
            bams_for_variant_calling.close()
            bams_for_variant_calling_tmp.close()
            bams_for_base_recalibration.close()
            bams_for_actual_recalibration.close()
        }
    })

// The next two lines are for vim to use the right syntax highlighting
// vim: syntax=groovy
// -*- mode: groovy;-*-
