#!/usr/bin/env nextflow

// This is a comment.
// vim: syntax=groovy
// -*- mode: groovy;-*-

// You should configure everything in the '*.nfconfig' file, not here.

input_fastqs = Channel
    .fromFilePairs(params.fastq_files_glob,size: -1)
    { file -> file.getBaseName() - ~/_n0[12]/ - ~/.fastq/ }
    .ifEmpty{ error "couldn't find those fastqs!" }

reference_prefix = Channel.fromPath(params.reference_prefix)

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
    echo true
    input:
        set val(pair_id), file(reads) from input_fastqs
        each reference_prefix
    output:
        set val(pair_id), file("aligned_reads.sam") into aligned_reads
    script:
    """
    module load ${params.modules.BWA}
    bwa mem -M -R "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.sequencing_platform}\\tPM:${params.sequencing_machine}\\tSM:${pair_id}" \
        ${reference_prefix} ${reads[0]} ${reads[1]} > \
        aligned_reads.sam
    """
}

process sort_and_bamize {
    echo true
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

aligned_reads_sorted.into { aligned_reads_sorted_metrics; aligned_reads_sorted_duplicates }

process collect_alignment_metrics {
    echo true
    publishDir "./output", mode: "copy"
    input:
        set val(pair_id), file(aligned_reads) \
            from aligned_reads_sorted_metrics
        each reference_prefix
    output:
        set val(pair_id), file("${pair_id}_alignment_metrics.txt"), \
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
        CollectAlignmentSummaryMetrics R=${reference_prefix} \
        I=${aligned_reads} O=${pair_id}_alignment_metrics.txt
    java -jar ${params.modules.PICARD_JAR} \
        CollectInsertSizeMetrics \
        INPUT=${aligned_reads} OUTPUT=${pair_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${pair_id}_insert_size_histogram.pdf 
    samtools depth -a ${aligned_reads} > ${pair_id}_depth_out.txt
    """
}

process mark_duplicates {
    echo true
    publishDir "./output", mode: "copy", pattern: "*_metrics_dedup.txt"
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
    echo true
    input:
        set val(pair_id), file(reads_dedup), file(reads_dedup_index) \
            from aligned_reads_deduplicated
        each reference_prefix
    output:
        set val(pair_id), val(1), \
            file("${pair_id}_realigned_reads.bam") \
            into realigned_around_indels
    script:
    """
    module load ${params.modules.GATK}
    java -jar ${params.modules.GATK_JAR} -T RealignerTargetCreator \
        -R ${reference_prefix} \
        -I ${reads_dedup} \
        -o ${pair_id}_realignment_targets.list
    java -jar ${params.modules.GATK_JAR} -T IndelRealigner \
        -R ${reference_prefix} \
        -I ${reads_dedup} \
        -targetIntervals ${pair_id}_realignment_targets.list \
        -o ${pair_id}_realigned_reads.bam"
    """
}

//process call_variants {
//    echo true
//    input:
//        set val(pair_id), val(round), file(input_bam)
//            from input_to_call_variants
//        each reference_prefix
//    output:
//        set val(pair_id), val(round), file("${pair_id}_variants.vcf")
//            into variants_called
//    script:
//    """
//    module load ${params.modules.GATK}
//    java -jar ${params.modules.GATK_JAR} -T HaplotypeCaller \
//        -R ${reference_prefix} \
//        -I ${input_bam} \
//        -o ${pair_id}_variants.vcf
//    """
//}

//	ROUND=$1
//	if [[ $ROUND -eq 1 ]];then
//		INPUT=${ID}_realigned_reads.bam
//		OUTPUT=${ID}_raw_variants.vcf
//		AFTEROK=$realignIndels
//	fi
//	if [[ $ROUND -eq 2 ]];then
//		INPUT=${ID}_recal_reads.bam
//		OUTPUT=${ID}_raw_variants_recal.vcf
//		AFTEROK=$applyBqsr
//	fi

//	
//	if [[ $ROUND -eq 1 ]];then
//		callVariants_1=$callVariants
//	fi
//	if [[ $ROUND -eq 2 ]];then
//		callVariants_2=$callVariants
//	fi
//}


//process extract_snps_and_indels {
//    echo true
//    input:
//        set val(pair_id), val(round), file(variants)
//            from variants_called
//        each reference_prefix
//    output:
//        set val(pair_id), val(round), file("${pair_id}_variants.vcf")
//            into variants_called
//    script:
//    """
//    module load ${params.modules.GATK}
//    java -jar ${params.modules.GATK_JAR} -T SelectVariants \
//        -R ${reference_prefix} \
//        -V ${variants} \
//        -selectType SNP \
//        -o ${pair_id}_snps_round${round}.vcf
//    java -jar ${params.modules.GATK_JAR} -T SelectVariants \
//        -R ${reference_prefix} \
//        -V ${variants} \
//        -selectType INDEL \
//        -o ${pair_id}_indel_round${round}.vcf
//    """
//}



//filter_snps(){
//        ROUND=$1
//        if [[ $ROUND -eq 1 ]];then
//                V=${ID}_raw_snps.vcf
//		O=${ID}_filtered_snps.vcf
//                AFTEROK=$extractSnps_1
//        fi
//        if [[ $ROUND -eq 2 ]];then
//                V=${ID}_raw_snps_recal.vcf
//                O=${ID}_filtered_snps_final.vcf
//                AFTEROK=$extractSnps_2
//        fi
//
//com="cd $FWD && \
//module load $GATK && \
//java -jar $GATK_JAR \
//-T VariantFiltration \
//-R $REF \
//-V $V \
//-filterName \"QD_filter\" \
//-filter \"QD < 2.0\" \
//-filterName \"FS_filter\" \
//-filter \"FS > 60.0\" \
//-filterName \"MQ_filter\" \
//-filter \"MQ < 40.0\" \
//-filterName \"SOR_filter\" \
//-filter \"SOR > 4.0\" \
//-o $O"
//response=\
//$(sbatch -J $ID.filterSnps$ROUND -o $ID.filterSnps$ROUND.out -e $ID.filterSnps$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
//stringarray=($response)
//filterSnps=${stringarray[-1]}
//echo $filterSnps >> $ID.log
//echo "FILTEREDSNPS: " $filterSnps
//echo "Filter SNPs Round $ROUND Submitted"
//
//        if [[ $ROUND -eq 1 ]];then
//                filterSnps_1=$filterSnps
//        fi
//        if [[ $ROUND -eq 2 ]];then
//                filterSnps_2=$filterSnps
//        fi
//}

//filter_indels(){
//        ROUND=$1
//        if [[ $ROUND -eq 1 ]];then
//                V=${ID}_raw_indels.vcf
//                O=${ID}_filtered_indels.vcf
//                AFTEROK=$extractSnps_1
//        fi
//        if [[ $ROUND -eq 2 ]];then
//                V=${ID}_raw_indels_recal.vcf
//                O=${ID}_filtered_indels_final.vcf
//                AFTEROK=$extractSnps_2
//        fi
//
//com="cd $FWD && \
//module load $GATK && \
//java -jar $GATK_JAR \
//-T VariantFiltration \
//-R $REF \
//-V $V \
//-filterName \"QD_filter\" \
//-filter \"QD < 2.0\" \
//-filterName \"FS_filter\" \
//-filter \"FS > 200.0\" \
//-filterName \"SOR_filter\" \
//-filter \"SOR > 10.0\" \
//-o $O"
//response=\
//$(sbatch -J $ID.filterIndels$ROUND -o $ID.filterIndels$ROUND.out -e $ID.filterIndels$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
//stringarray=($response)
//filterIndels=${stringarray[-1]}
//echo $filterIndels >> $ID.log
//echo "FILTEREDINDELS: " $filterIndels
//echo "Filter Indels Round $ROUND Submitted"
//
//}

//do_bqsr(){
//	#todo: knownSites input shouldnt be full raw_variants.vcf file but only the TOP variants!
//	ROUND=$1
//	if [[ $ROUND -eq 1 ]];then
//		POST=''
//		OUT=${ID}_recal_data.table
//		AFTEROK="$filterSnps_1:$filterIndels"
//	fi
//	if [[ $ROUND -eq 2 ]];then
//		POST='-BQSR '${ID}_recal_data.table
//		OUT=${ID}_post_recal_data.table
//		AFTEROK=$bqsr_1
//	fi
//
//com="cd $FWD && \
//module load $GATK && \
//java -jar $GATK_JAR \
//-T BaseRecalibrator \
//-R $REF \
//-I ${ID}_realigned_reads.bam \
//-knownSites ${ID}_filtered_snps.vcf \
//-knownSites ${ID}_filtered_indels.vcf \
//$POST \
//-o $OUT"
//response=\
//$(sbatch -J $ID.bqsr$ROUND -o $ID.bqsr$ROUND.out -e $ID.bqsr$ROUND.err --dependency=afterok:$AFTEROK --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
//stringarray=($response)
//bqsr=${stringarray[-1]}
//echo $bqsr >> $ID.log
//echo "BQSR: " $bqsr
//echo "BQSR Round $ROUND Submitted"
//
//	if [[ $ROUND -eq 1 ]];then
//		bqsr_1=$bqsr
//	fi
//	if [[ $ROUND -eq 2 ]];then
//		bqsr_2=$bqsr
//	fi
//}

//analyze_covariates(){
//com="cd $FWD && \
//module load $GATK && \
//module load $R && \
//java -jar $GATK_JAR \
//-T AnalyzeCovariates \
//-R $REF \
//-before ${ID}_recal_data.table \
//-after ${ID}_post_recal_data.table \
//-plots ${ID}_recalibration_plots.pdf"
//response=\
//$(sbatch -J $ID.analyzeCovariates -o $ID.analyzeCovariates.out -e $ID.analyzeCovariates.err --dependency=afterok:$bqsr_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
//stringarray=($response)
//analyzeCovariates=${stringarray[-1]}
//echo $analyzeCovariates >> $ID.log
//echo "ANALYZECOVARIATES: " $analyzeCovariates
//echo "AnalyzeCovariates Submitted"
//}

//apply_bqsr(){
//com="cd $FWD && \
//rm ${ID}_dedup_reads.bam ${ID}_dedup_reads.bai && \
//module load $GATK && \
//java -jar $GATK_JAR \
//-T PrintReads \
//-R $REF \
//-I ${ID}_realigned_reads.bam \
//-BQSR ${ID}_recal_data.table \
//-o ${ID}_recal_reads.bam"
//response=\
//$(sbatch -J $ID.applyBqsr -o $ID.applyBqsr.out -e $ID.applyBqsr.err --dependency=afterok:$bqsr_2 --kill-on-invalid-dep=yes --mail-user=$EMAIL --mail-type=FAIL --nodes=1 -t 4:00:00 --mem=60000 --wrap="$com")
//stringarray=($response)
//applyBqsr=${stringarray[-1]}
//echo $applyBqsr >> $ID.log
//echo "APPLYBQSR: " $applyBqsr
//echo "Apply BQSR Submitted"
//}

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
