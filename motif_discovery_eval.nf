#!/usr/bin/env nextflow

// params.inputDir1 = "${HOME}/nextflow_test"
// params.inputDir2 = "${HOME}/nextflow_test"
// params.ywuv_htselex = "/home/hughespub/SELEX_Data/HT-SELEX_Fastqs_Ready_For_SRA_etc/YWUV_and_RoziHadiAttackATAC_240711_A00546_0178_AH57KMDRX5"
// params.ywt_htselex = "/home/hughespub/finishTF/YWT_Mostly_Isaac/FinalFastq_Trimmed_and_Filtered_Ready"
// params.yww_htselex = "/home/hughespub/SELEX_Data/HT-SELEX_Fastqs_Ready_For_SRA_etc/HTSELEX_BcorrectFilteredTrimmed"
// // params.ght_input1 = "/home/hughespub/SELEX_Data/GHT-SELEX_Fastqs_Ready_For_SRA_etc/YWUV_and_RoziHadiAttackATAC_240711_A00546_0178_AH57KMDRX5"
// params.ght_input = "/home/hughespub/SELEX_Data/GHT-SELEX_Fastqs_Ready_For_SRA_etc"
// params.ght_peaks = "/home/hughespub/ahcorcha/Transposone_TFs_MAGIX_peaks_08_11_24"
// params.outputDir = "${HOME}/TEHMM_proj/analysis/selex_results/motif_pipeline"

include { MOTIF_BENCHMARK; BENCHMARK_PREP; LOGO_MAKER; BENCHMARK_PREP_PEAKS; MOTIF_BENCHMARK_PEAKS } from './motif_benchmark.nf'
include { MOTIF_DISCOVERY; RUN_BEESEM; RUN_DIMONT_HTS; RUN_DIMONT_GHT } from './motif_finders.nf'
include { JOIN_GHT_READS; RUN_BOWTIE; PEAK_CALL; KNEEDLE_PEAKS; SPLIT_PEAKS; EXTRACT_TOP_PEAKS } from './peak_processing.nf'


process CAT_FASTQ {
    conda "${HOME}/micromamba"
    input:
        tuple val(basenames), val(cycles), val(exp_id), path(fastqs) // to deal with ght fastqs with identical names since they're the read pair joining step

    // output the merged fastq file, and the exp_id with "merged" appended
    output:
        tuple path("merged.fastq.gz"), val("${exp_id}_merged")
    script:
        """
        gunzip -c ${fastqs[0]} ${fastqs[1]} ${fastqs[2]} > merged.fastq
        gzip merged.fastq
        """
}

process DIMONT_HT_PREP {
    conda "${HOME}/micromamba"
    input:
        tuple val(basenames), val(cycles), val(exp_id), path(fastqs)
    
    // output the merged fastq file, and the exp_id with "merged_dimont" appended
    output:
        tuple path("merged_dimont.fa.gz"), val("${exp_id}_merged_dimont")

    script:
        """
        ## extract equal number of random reads from each (190k to account for 
        ## subsequent split into test and train)
        dimont_hts_format.sh 1 ${fastqs[0]} > merged_dimont.fa
        dimont_hts_format.sh 2 ${fastqs[1]} >> merged_dimont.fa
        dimont_hts_format.sh 3 ${fastqs[2]} >> merged_dimont.fa
        gzip merged_dimont.fa
        """
}

// Define the process to deduplicate the fastq files
process DE_DUP {
    label 'dedupe'
    conda "${HOME}/micromamba"
    input:
        tuple path(fastq), val(baseName)

    output:
        tuple path("dedupe.fastq.gz"), val(baseName)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        # dedupe.sh in=$fastq out="dedupe.fastq.gz"
        java -ea -Xmx10939m -Xms10939m \
        -cp /home/hugheslab1/iyellan/micromamba/opt/bbmap-39.08-0/current/ \
        jgi.Dedupe in=$fastq out=dedupe.fastq.gz
        """
}

// Use seqkit to split into 70% train and 30% test
process SPLIT {
    conda "${HOME}/micromamba"
    label 'high_mem'
    input:
        tuple path(fastq), val(baseName)

    output:
        tuple path("train.fastq.gz"), val(baseName), emit: train
        tuple path("test.fastq.gz"), val(baseName), emit: test

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true
    
    script:
        """
        seqkit shuffle --rand-seed 1 $fastq -o shuffled.fastq.gz
        seqkit sample -p 0.3 shuffled.fastq.gz -o test.fastq.gz
        seqkit seq -n -i test.fastq.gz -o test_ids.txt
        seqkit grep -v -f test_ids.txt shuffled.fastq.gz -o train.fastq.gz
        """
}

process SPLIT_FASTA {
    conda "${HOME}/micromamba"
    label 'high_mem'
    input:
        tuple path(fasta), val(baseName)

    output:
        tuple path("train.fasta.gz"), val(baseName), emit: train
        tuple path("test.fasta"), val(baseName), emit: test

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        gunzip -c $fasta | $HOME/scripts/fasta_linearize.sh \
        > shuffed_fasta.txt
        n_seqs=\$(cat shuffed_fasta.txt | wc -l)

        n_train=\$(printf "%.0f" \$(echo "\$n_seqs * 0.7" | bc -l))
        n_test=\$((n_seqs-n_train))

        head -n \$n_train shuffed_fasta.txt | awk 'BEGIN{OFS="\\n"} {print ">"\$1,\$2}' \
        | sed 's/,/ /g' > train.fasta
        gzip -f train.fasta
        tail -n \$n_test shuffed_fasta.txt | awk 'BEGIN{OFS="\\n"} {print ">"\$1,\$2}' \
        | sed 's/,/ /g' > test.fasta # for some reason pwm-eval doesn't work if the dimont fasta is gzipped
        """
}


// Use the GRECO-BIT folk's k-mer enrichment script to subset the training data fastq 
// to the 10,000 sequences with the highest k-mer enrichment scores at both 5 and 10 kmer lengths
process KMER_ENRICH {
    label 'kmer_enr'
    conda "${HOME}/micromamba"
    input:
        tuple path(fastq), val(baseName), val(kmerLen)

    output:
        tuple path("kmer_enr_${kmerLen}.fasta.gz"), val(baseName), val(kmerLen)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        $HOME/software/HT-SELEX-kmer-filtering/extract_topk.py -d -k ${kmerLen} \
        --fastq $fastq --out kmer_enr_${kmerLen}.fasta.gz -b 10000 -a 0
        """
}

process KMER_ENRICH_MULTITHREAD {
    label 'kmer_enr_multithread'
    conda '/home/hugheslab1/iyellan/micromamba'
    input:
        tuple path(fastq), val(baseName), val(kmerLen)

    output:
        tuple path("kmer_enr_${kmerLen}.fasta.gz"), val(baseName), val(kmerLen)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        ## account for ght fragments having different lengths (-m flag)
        $HOME/software/HT-SELEX-kmer-filtering/extract_topk_multithread.py \
        -m -d -k ${kmerLen} --fastq $fastq --out kmer_enr_${kmerLen}.fasta.gz \
        -b 10000 -a 0
        """
}

process SEQ_PREP {
    conda '/home/hugheslab1/iyellan/micromamba'
    label 'seq_prep'
    input:
        tuple path(fasta), val(baseName), val(kmerLen)

    output:
        // tuple path("unzipped.fasta"), val(baseName), emit: unzipped
        tuple path("unzipped_${kmerLen}_unique_ids.fasta"), val(baseName), val(kmerLen)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        gunzip -c "$fasta" | sed 's/ //g' > unzipped_"${kmerLen}".fasta
        cat unzipped_"${kmerLen}".fasta \
        | fasta-unique-names > unzipped_"${kmerLen}"_unique_ids.fasta
        """

}

// run DE_DUP and SPLIT as workflow
workflow HTSELEX {
    take:
        inputDir
    
    main:
        ////   HT READ PROCESSING /////////
        // Channel for individual HT-SELEX files
        Channel
        .fromPath( ["${params.ywuv_htselex}/*.fastq.gz", "${params.ywt_htselex}/*40N*_pTH*.fastq.gz", 
        "${params.ywt_htselex}/*40N*_UT380*.fastq.gz", "${params.yww_htselex}/*.fastq.gz"] )
        .map {fastq ->
            def baseName = fastq.getSimpleName()
            tuple(fastq, baseName)
        }
        .set { ht_fastq_tuples_individual }

        // Channel for YWU/V HT-SELEX fastq files that are to be catted
        ywuv_for_catting = Channel
            .fromPath( ["${params.ywuv_htselex}/*.fastq.gz"] )
            .map { file -> 
                // Extract cycle (A_1, A_2, etc.) and unique identifier (e.g., pTH14647_TT40NCTCGTC_eGFP_IVT_SXXX)
                def baseName = file.getSimpleName()
                def cycle = baseName.split('_')[1..2].join('_')  // e.g., A_1
                def exp_id = baseName.split('_')[0,3..5].join('_') // excluding the cycle and SXXX

                tuple(baseName, cycle, exp_id, file)
            }

        // Channel for YWT fastq files to be catted; requires different processing
        ywt_for_catting = Channel
            .fromPath( ["${params.ywt_htselex}/*40N*_pTH*.fastq.gz", 
            "${params.ywt_htselex}/*40N*_UT380*.fastq.gz"] )
            .map { file -> 
                // Extract cycle (A_1, A_2, etc.) and unique identifier (e.g., pTH14647_TT40NCTCGTC_eGFP_IVT_SXXX)
                def baseName = file.getSimpleName()
                def cycle = baseName.split('_')[4]  // e.g., A_1
                def exp_id = baseName.split('_')[0..3].join('_') // excluding the cycle and SXXX

                tuple(baseName, cycle, exp_id, file)
            }
        
        // Channel for YWW HT-SELEX fastq files to be catted; requires different processing
        yww_for_catting = Channel
            .fromPath( ["${params.yww_htselex}/*_pTH*.fastq.gz"] )
            .map { file -> 
                // Extract cycle (A_1, A_2, etc.) and unique identifier (e.g., pTH14647_TT40NCTCGTC_eGFP_IVT_SXXX)
                def baseName = file.getSimpleName()
                def cycle = baseName.split('_')[1]  // e.g., A_1
                def exp_id = baseName.split('_')[2..5].join('_') // excluding the cycle and SXXX

                tuple(baseName, cycle, exp_id, file)
            }
        // Combine the two HT-SELEX catting channels
        ht_selex_for_catting = ywuv_for_catting.mix(ywt_for_catting).mix(yww_for_catting)

        // combine the merged GHT fastq files with the ht-selex fastq files for catting
        ht_ght_selex_for_catting = ht_selex_for_catting

        // also create separate ght channel of individual ght fastqs for deduplication
        // ght_individual = ght_merged.map { baseName, cycle, exp_id, merged_fq ->
        //     tuple(merged_fq, baseName)
        // }
        //// Create catted fastq files for HT-SELEX and GHT-SELEX across cycles ////
        // group by exp_id
        ht_ght_selex_for_catting_grouped = ht_ght_selex_for_catting.groupTuple(by:2)
        // concatenate the fastqs for each group
        catted_fqs = CAT_FASTQ(ht_ght_selex_for_catting_grouped)
        
        // also concatenate ht selex reads separately, for motif finding with dimont-hts
        ht_selex_for_catting_grouped = ht_selex_for_catting.groupTuple(by:2)
        dimont_hts_fas = DIMONT_HT_PREP(ht_selex_for_catting_grouped)

        //// Handling individual cycles
        // deduplicate the HT-SELEX fastqs for individual cycles
        deduped_ht_selex = DE_DUP(ht_fastq_tuples_individual)

        // mix the deduped fastq_tuples with the ght_individuals into a single channel
        ht_ght_individual = deduped_ht_selex//.mix(ght_individual)

        // finally, mix the individual cycle fqs with the catted_fqs into a single channel,
        // because the concatenated fastq files are expected to have duplicate sequences across cycles
        // so bypassing deduplication for the concatenated fastqs and GHT fastqs
        deduped_and_merged = ht_ght_individual.mix(catted_fqs)

        // split the fastqs into training and testing sets
        split = SPLIT(deduped_and_merged)

        // also split the dimont fastas
        dimont_hts_split = SPLIT_FASTA(dimont_hts_fas)

        // find the 10,000 sequences with the highest k-mer enrichment scores in the training set
        //// create a kmer length channel
        kmer_lens = Channel.of(5, 8)
        train_kmer_lens = split.train.combine(kmer_lens)

        // divide train_kmer_lens channel based on if basename contains "GHT" or not
        ht_train_kmer_lens = train_kmer_lens.filter { fastq, baseName, kmerLen ->
            !baseName.contains("GHT")
        }
        // run the kmer enrichment process
        kmerEnr = KMER_ENRICH(ht_train_kmer_lens)
        unzipped = SEQ_PREP(kmerEnr)

        motifs = MOTIF_FINDERS(unzipped)
        // dimont needs to be run separately, because it uses all reads as input
        dimont = MOTIFS_DIMONT(dimont_hts_split)
        // mix the motif finding results together
        motif_dir_tuples = motifs.mix(dimont)

        // benchmark the motifs
        MOTIF_BENCHMARK_WORKFLOW(motif_dir_tuples)
}

workflow GHTPEAKCALLING {
    take:
        inputDir
    
    main:
        println "Input directory: ${inputDir}"

        Channel
            .fromFilePairs("${inputDir}/GHT0*_Cyc*Human_*_YWW_*R{1,2}_001.fastq.gz", flat: true)
            .map { grouping_key, pair1, pair2  -> 
                // Extract cycle (A_1, A_2, etc.) and unique identifier (e.g., pTH14647_TT40NCTCGTC_eGFP_IVT_SXXX)
                def baseName = pair1.getSimpleName()
                def cycle = baseName.split('_')[1].join('_')  // e.g., A_1
                def exp_id = baseName.split('_')[0,2..6].join('_') // excluding the cycle and SXXX
                tuple(baseName, pair1, pair2)
            }
            .set { ght_fastq_tuples }
        // CHANNEL FOR GHT CONTROLS
        Channel
        .fromFilePairs("${inputDir}/*Naive*Human*R{1,2}.fastq.gz", flat: true)
        .map { grouping_key, pair1, pair2  -> 
            // Extract cycle (A_1, A_2, etc.) and unique identifier (e.g., pTH14647_TT40NCTCGTC_eGFP_IVT_SXXX)
            def baseName = pair1.getSimpleName()
            tuple("GHT_HUMAN_CONTROL", pair1, pair2)
        }
        .set { ght_control_fastq_tuples }
        
        // MAP GHT FASTQ FILES TO THE HUMAN GENOME WITH BOWTIE
        ght_mapped_bams = RUN_BOWTIE(ght_fastq_tuples.mix(ght_control_fastq_tuples))

        // SEPARATE CONTROLS FROM EXPERIMENTS
        ght_mapped_bam = ght_mapped_bams.filter { bam, baseName ->
            !baseName.contains("GHT_HUMAN_CONTROL")
        }
        ght_control_bam = ght_mapped_bams.filter { bam, baseName ->
        baseName.contains("GHT_HUMAN_CONTROL")
        }

        // COMBINE ALL BAM FILES WITH THE SAME BASENAME INTO A TUPLE
        ght_mapped_bam_tuples = ght_mapped_bam.groupTuple(by:1)
        ght_control_bam_tuples = ght_control_bam.groupTuple(by:1)

        // CALL PEAKS WITH MACS3
        ght_macse_peaks = PEAK_CALL(ght_mapped_bam_tuples.combine(ght_control_bam_tuples))

    emit:
        ght_macse_peaks
}

workflow GHTSELEX {
    ////   GHT READ PROCESSING - CALL PEAKS, SPLIT PREAKS, EXTRACT PEAK SEQS /////////
    // Create channel of GHT fastq files, associate pairs of files with the same exp_id
    take:
        inputDir
    main:
        // if input is ght fastq, run peak calling
        if("${params.inputType}" == "fastq"){
            ght_peak_beds = GHTPEAKCALLING(inputDir)
        }
        else if("${params.inputType}" == "peaks"){
            Channel
            .fromPath("${params.inputDir}/*/*_LTR_results_all_with_eFDR.bed")
            .map { bed ->
            def baseName = bed.getSimpleName() + "_PEAKS"
            def peakSetType = "full"
            tuple(bed, baseName, peakSetType)
            }
            .set { full_peaks }

            // peak file processing
            // use kneedle method to subset the ght peak bed files to "real" peaks
            kneedle_peaks = KNEEDLE_PEAKS(full_peaks)

            // mix the kneedle and full peak bed files
            ght_peak_beds = kneedle_peaks.mix(full_peaks)
        }
        
        // Run the JOIN_GHT_READS process on the GHT fastq files
        // ght_merged = JOIN_GHT_READS(ght_fastq_tuples)
        
        // split the ght peak bed file into training and testing sets
        // do for both kneedle-filtered (for most motif finders) and full peak sets (for dimont)
        // macse called peaks can be run on all the motif finders
        split_peaks_all = SPLIT_PEAKS(ght_peak_beds)

        split_peaks_train = split_peaks_all.train_bed.filter { bed, baseName, peakSetType ->
            peakSetType =="kneedle" || peakSetType == "macs"
        }
        split_peaks_test = split_peaks_all.test_bed.filter { bed, baseName, peakSetType ->
            peakSetType =="kneedle" || peakSetType == "macs"
        }
        dimont_split_peaks_train = split_peaks_all.train_bed.filter { bed, baseName, peakSetType ->
            peakSetType =="full" || peakSetType == "macs"
        }
        dimont_split_peaks_test = split_peaks_all.test_bed.filter { bed, baseName, peakSetType ->
            peakSetType =="full" || peakSetType == "macs"
        }

        // grab top scoring peaks from the ght peak bed file, extract fastas
        n_peaks = Channel.of(50, 200, 500)
        ght_peaks_topn = split_peaks_train.combine(n_peaks)
        ght_peaks_topn_fas = EXTRACT_TOP_PEAKS(ght_peaks_topn)

        motifs = MOTIF_FINDERS(ght_peaks_topn_fas)
        // run dimont ght on the ght peaks (also extracts fasta sequences)
        dimont = MOTIFS_DIMONT(dimont_split_peaks_train)
        motif_dir_tuples = motifs.mix(dimont)

        // benchmark the motifs
        test_set = split_peaks_test.mix(dimont_split_peaks_test)
        MOTIF_BENCHMARK_WORKFLOW(motif_dir_tuples, test_set)
}

workflow MOTIFS_DIMONT {
    take: dimontIn
    
    main:
        if(params.inputExp == "HT") {
            dimont = RUN_DIMONT_HTS(dimontIn)
        }
        else if(params.inputExp == "GHT") {
            dimont = RUN_DIMONT_GHT(dimontIn)
        }
    emit: 
        dimont
}

workflow MOTIF_FINDERS {
    take: fastas

    main:
        // create channel of motif caller types, combine with unzipped channel as 5th item in tuple
        motif_callers = Channel.of( "meme", "streme", "homer" )
        fastas_motif_callers = fastas.combine(motif_callers)
        motifs = MOTIF_DISCOVERY(fastas_motif_callers)
        // beesem needs to be run separately because it requires a different conda environment
        beesem = RUN_BEESEM(fastas_motif_callers)

        // mix the motif finding results together
        motif_dir_tuples = motifs.mix(beesem)

    emit:
        motif_dir_tuples
}

workflow MOTIF_BENCHMARK_WORKFLOW {
    take: 
        motif_dir_tuples
        test_set

    main:
        // prepare the benchmarking process
        //// channel of fractions of positives to use for evaluation
        top_fracs = Channel.of(0.01, 0.1, 0.5)

        pfm_tuples = motif_dir_tuples.map { pfm_dir, base, kmerLen ->
            // Find the .pfm files in the MEME output directory
            def fils = file("${pfm_dir}/*.pfm").collect()
            // Pair each .pfm file with the base name
            fils.collect { fil -> tuple(base, fil, kmerLen) }
        }.flatten().collate(3)

        if("${params.inputExp}" == "HT") {
            //// create positives and negatives (dinuc shuffle) from read test set
            prep = BENCHMARK_PREP(test_set)
        }
        else if("${params.inputExp}" == "GHT") {
            // top fracs are needed earlier in the peak prep stage whereas they're used in the read evaluation stage
            prep = BENCHMARK_PREP_PEAKS(test_set.combine(top_fracs))
        }

        // prep = BENCHMARK_PREP(split.test.mix(dimont_hts_split.test))
        // prep_peaks = BENCHMARK_PREP_PEAKS(split_peaks_test.mix(dimont_split_peaks_test).combine(top_fracs))
        

        // filter pfm_tuples to read-based and peak-based motifs
        // pfm_tuples_reads = pfm_tuples.filter { base, pfm, kmerLen ->
        //     !base.contains("PEAKS")
        // }
        // pfm_tuples_peaks = pfm_tuples.filter { base, pfm, kmerLen ->
        //     base =~ /PEAKS/
        // }

        // Combine motif results with prep output (pairing based on order of emission)
        pfm_bmark_prep_tuples = pfm_tuples.combine(prep, by: 0)
        // pfm_bmark_prep_tuples_peaks = pfm_tuples_peaks.combine(prep_peaks, by: 0)
        // combine pfm_bmark_prep_top_frac_tuples_reads with a channel of top fractions (already done for peaks)

        if("${params.inputExp}" == "HT") {
            MOTIF_BENCHMARK(pfm_bmark_prep_tuples)

        }
        else if("${params.inputExp}" == "GHT") {
            MOTIF_BENCHMARK_PEAKS(pfm_bmark_prep_tuples)

        }
        // Run the benchmarking process with the combined motifs
        // create logos from the pfm files
        LOGO_MAKER(pfm_tuples)
}

workflow {
    /// MAIN WORFKLOW ///
    if("${params.inputExp}" == "HT") {
        HTSELEX("${params.inputDir}")
    }
    else if("${params.inputExp}" == "GHT") {
        GHTSELEX("${params.inputDir}")
    }
}
