#!/usr/bin/env nextflow

// params.inputDir1 = "${HOME}/nextflow_test"
// params.inputDir2 = "${HOME}/nextflow_test"
params.ywuv_htselex = "/home/hughespub/SELEX_Data/HT-SELEX_Fastqs_Ready_For_SRA_etc/YWUV_and_RoziHadiAttackATAC_240711_A00546_0178_AH57KMDRX5"
params.ywt_htselex = "/home/hughespub/finishTF/YWT_Mostly_Isaac/FinalFastq_Trimmed_and_Filtered_Ready"
params.ght_input = "/home/hughespub/SELEX_Data/GHT-SELEX_Fastqs_Ready_For_SRA_etc/YWUV_and_RoziHadiAttackATAC_240711_A00546_0178_AH57KMDRX5"
params.ght_peaks = "/home/hughespub/ahcorcha/Transposone_TFs_MAGIX_peaks_08_11_24"
params.outputDir = "${HOME}/TEHMM_proj/analysis/selex_results/motif_pipeline"

include { MOTIF_BENCHMARK; BENCHMARK_PREP; LOGO_MAKER; BENCHMARK_PREP_PEAKS; MOTIF_BENCHMARK_PEAKS } from './motif_benchmark.nf'

process JOIN_GHT_READS {
    conda '/home/hugheslab1/iyellan/micromamba'
    label 'merge_reads'
    input:
        tuple val(baseName), val(cycle), val(exp_id), path(pair1), path(pair2)
    output:
        tuple val(baseName), val(cycle), val(exp_id), path("ght_merged_${cycle}.fastq.gz")
    script:
        """
        bbmerge-auto.sh in1=$pair1 in2=$pair2 out=ght_merged_${cycle}.fastq.gz rem k=62 extend2=50 ecct \
        -Xmx10939m
        """
}

process CAT_FASTQ {
    conda '/home/hugheslab1/iyellan/micromamba'
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

// Define the process to deduplicate the fastq files
process DE_DUP {
    label 'dedupe'
    conda '/home/hugheslab1/iyellan/micromamba'
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

// Use reformat.sh to convert fastq to fasta and split into 70% train and 30% test
process SPLIT {
    conda '/home/hugheslab1/iyellan/micromamba'
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

process SPLIT_PEAKS {
    conda '/home/hugheslab1/iyellan/micromamba'
    input:
        tuple path(bed), val(baseName)

    output:
        tuple path("train.bed"), val(baseName), emit: train_bed
        tuple path("test.bed"), val(baseName), emit: test_bed
    
    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        n_lines=\$(cat $bed | wc -l)
        pos_lines=\$(printf "%.0f" \$(echo "\$n_lines * 0.7" | bc -l))
        # round up to nearest integer
        neg_lines=\$((n_lines - pos_lines)) 
        shuf <(tail -n +2 $bed) > shuffled.bed
        head -n \$pos_lines shuffled.bed > train.bed
        tail -n \$neg_lines shuffled.bed > test.bed
        """
}

process EXTRACT_TOP_PEAKS {
    conda '/home/hugheslab1/iyellan/micromamba'
    //nPeaks in the output can serve as a placeholder instead of kmerLen
    input:
        tuple path(shuff_bed), val(baseName), val(nPeaks)

    output:
        tuple path("top_${nPeaks}.fasta"), val(baseName), val(nPeaks)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        sort -k5 -gr $shuff_bed | head -n $nPeaks \
        | awk 'BEGIN{OFS=FS="\t"} {print \$1,\$2,\$3,\$4,\$5}' > top_${nPeaks}.bed
        bedtools getfasta -fi \$HOME/data/hg38.fa \
        -bed top_${nPeaks}.bed -fo top_${nPeaks}.fasta
        """
}

// Use the GRECO-BIT folk's k-mer enrichment script to subset the training data fastq 
// to the 10,000 sequences with the highest k-mer enrichment scores at both 5 and 10 kmer lengths
process KMER_ENRICH {
    label 'kmer_enr'
    conda '/home/hugheslab1/iyellan/micromamba'
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

process MOTIF_DISCOVERY {
    label 'multithread'
    conda '/home/hugheslab1/iyellan/micromamba'

    input:
        tuple path(fasta), val(baseName), val(kmerLen), val(motifCaller)

    output:
        tuple path("${motifCaller}_out_${kmerLen}/"), val(baseName), val(kmerLen)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        if("${motifCaller}" == "meme")
            """
            meme $fasta -oc meme_out_${kmerLen}/ -minw 5 -maxw 20 -p 4 -dna -nmotifs 3
            $HOME/TEHMM_proj/selex_motif_scripts/meme_xml_parse.py meme_out_${kmerLen}/meme.xml
            """
        
        else if("${motifCaller}" == "streme")
            """
            streme --p $fasta --oc streme_out_${kmerLen}/ --minw 5 --maxw 20 --neval 40 --nref 8 --niter 25 --dna
            $HOME/TEHMM_proj/selex_motif_scripts/meme_xml_parse.py streme_out_${kmerLen}/streme.xml
            """
        
        else if("${motifCaller}" == "homer")
            """
            fasta-shuffle-letters -kmer 2 -seed 42 $fasta > dummy_file
            homer2 denovo -i $fasta -p 4 -b dummy_file -len 6 -o homer_out_6.txt
            homer2 denovo -i $fasta -p 4 -b dummy_file -len 10 -o homer_out_10.txt
            homer2 denovo -i $fasta -p 4 -b dummy_file -len 14 -o homer_out_14.txt
            for i in 6 10 14; do
                n_lines=\$((i+1))
                head -n \${n_lines} homer_out_\${i}.txt | tail -n +2 > homer_out_\${i}_1.pfm
                head -n \$((n_lines*2)) homer_out_\${i}.txt | tail -n \$((n_lines - 1)) > homer_out_\${i}_2.pfm
            done
            mkdir -p homer_out_${kmerLen}
            mv homer_out_*.txt homer_out_*.pfm homer_out_${kmerLen}/
            """
}

process RUN_BEESEM {
    conda '/home/hugheslab1/iyellan/micromamba/envs/py27_env'
    label 'RUN_BEESEM'
    input:
        tuple path(fasta), val(baseName), val(kmerLen)
    
    output:
        tuple path("beesem_out_${kmerLen}/"), val(baseName), val(kmerLen), emit: beesem
    
    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        $HOME/scripts/fasta_linearize.sh $fasta \
        | awk -F"\t" '{print \$2}' | uniq -c | awk 'BEGIN{OFS="\t"} {print toupper(\$2),\$1}' \
        > beesem_formatted.txt
        python $HOME/software/BEESEM/beesem.py -o beesem_out_${kmerLen} "${baseName}" beesem_formatted.txt
        ## parse the output
        beesem_pfm=\$(ls beesem_out_${kmerLen}/"${baseName}"_rep=1_phs=10/results/pfm_*.txt)
        tail -n +3 \${beesem_pfm} > beesem_out_${kmerLen}/beesem_rfmt.pfm
        """
}



// run DE_DUP and SPLIT as workflow
workflow {
     ////   HT READ PROCESSING /////////
    // Channel for individual HT-SELEX files
    Channel
    .fromPath( ["${params.ywuv_htselex}/*.fastq.gz", "${params.ywt_htselex}/*40N*_pTH*.fastq.gz", 
    "${params.ywt_htselex}/*40N*_UT380*.fastq.gz"] )
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
    // Combine the two HT-SELEX catting channels
    ht_selex_for_catting = ywuv_for_catting.mix(ywt_for_catting)

    ////   GHT READ PROCESSING /////////
    // Create channel of GHT fastq files, associate pairs of files with the same exp_id
    Channel
    .fromFilePairs("${params.ght_input}/YW*GHT_Human_eGFP_IVT_*_R{1,2}_001.fastq.gz", flat: true)
    .map { grouping_key, pair1, pair2  -> 
        // Extract cycle (A_1, A_2, etc.) and unique identifier (e.g., pTH14647_TT40NCTCGTC_eGFP_IVT_SXXX)
        def baseName = pair1.getSimpleName()
        def cycle = baseName.split('_')[1..2].join('_')  // e.g., A_1
        def exp_id = baseName.split('_')[0,3..5].join('_') // excluding the cycle and SXXX

        tuple(baseName, cycle, exp_id, pair1, pair2)
    }
    .set { ght_fastq_tuples }
    
    // Run the JOIN_GHT_READS process on the GHT fastq files
    ght_merged = JOIN_GHT_READS(ght_fastq_tuples)
    
    // GHT peak file processing
    Channel
    .fromPath("${params.ght_peaks}/*/*_LTR_results_all_with_eFDR.bed")
    .map { bed ->
        def baseName = bed.getSimpleName() + "_PEAKS"
        tuple(bed, baseName)
    }
    .set { ght_peak_beds }

    // combine the merged GHT fastq files with the ht-selex fastq files for catting
    ht_ght_selex_for_catting = ht_selex_for_catting.mix(ght_merged)

    // also create separate ght channel of individual ght fastqs for deduplication
    ght_individual = ght_merged.map { baseName, cycle, exp_id, merged_fq ->
        tuple(merged_fq, baseName)
    }

    //// Create catted fastq files for HT-SELEX and GHT-SELEX across cycles ////
    // group by exp_id
    ht_ght_selex_for_catting_grouped = ht_ght_selex_for_catting.groupTuple(by:2)
    // concatenate the fastqs for each group
    catted_fqs = CAT_FASTQ(ht_ght_selex_for_catting_grouped)

    //// Handling individual cycles
    // deduplicate the HT-SELEX fastqs for individual cycles
    deduped_ht_selex = DE_DUP(ht_fastq_tuples_individual)

    // mix the deduped fastq_tuples with the ght_individuals into a single channel
    ht_ght_individual = deduped_ht_selex.mix(ght_individual)

    // finally, mix the individual cycle fqs with the catted_fqs into a single channel,
    // because the concatenated fastq files are expected to have duplicate sequences across cycles
    // so bypassing deduplication for the concatenated fastqs and GHT fastqs
    deduped_and_merged = ht_ght_individual.mix(catted_fqs)

    // split the fastqs into training and testing sets
    split = SPLIT(deduped_and_merged)

    // split the ght peak bed file into training and testing sets
    split_peaks = SPLIT_PEAKS(ght_peak_beds)
    // find the 10,000 sequences with the highest k-mer enrichment scores in the training set
    //// create a kmer length channel
    kmer_lens = Channel.of(5, 8)
    train_kmer_lens = split.train.combine(kmer_lens)

    // divide train_kmer_lens channel based on if basename contains "GHT" or not
    ght_train_kmer_lens = train_kmer_lens.filter { fastq, baseName, kmerLen ->
        baseName =~ /GHT/
    }
    ht_train_kmer_lens = train_kmer_lens.filter { fastq, baseName, kmerLen ->
        !baseName.contains("GHT")
    }
    // run the kmer enrichment process
    kmerEnr_ht = KMER_ENRICH(ht_train_kmer_lens)
    kmerEnr_ght = KMER_ENRICH_MULTITHREAD(ght_train_kmer_lens)
    kmerEnr = kmerEnr_ht.mix(kmerEnr_ght)
    unzipped = SEQ_PREP(kmerEnr)

    // grab top scoring peaks from the ght peak bed file, extract fastas
    n_peaks = Channel.of(50, 200, 500)
    ght_peaks_topn = split_peaks.train_bed.combine(n_peaks)
    ght_peaks_topn_fas = EXTRACT_TOP_PEAKS(ght_peaks_topn)

    // mix the peak fasta files with the read fastas
    unzipped_ght_peaks = unzipped.mix(ght_peaks_topn_fas)

    // create channel of motif caller types, combine with unzipped channel as 5th item in tuple
    motif_callers = Channel.of( "meme", "streme", "homer" )
    unzipped_motif_callers = unzipped_ght_peaks.combine(motif_callers)
    motifs = MOTIF_DISCOVERY(unzipped_motif_callers)
    // beesem needs to be run separately because it requires a different conda environment
    beesem = RUN_BEESEM(unzipped_ght_peaks)

    // prepare the benchmarking process
    top_fracs = Channel.of(0.01, 0.1, 0.5)
    prep = BENCHMARK_PREP(split.test)
    prep_peaks = BENCHMARK_PREP_PEAKS(split_peaks.test_bed.combine(top_fracs)) // top fracs are needed earlier in the peak prep stage
    // mix the .pfm files from MEME, STREME, HOMER, and BEESEM into a single channel
    motif_dir_tuples = motifs.mix(beesem)
    pfm_tuples = motif_dir_tuples.map { pfm_dir, base, kmerLen ->
        // Find the .pfm files in the MEME output directory
        def fils = file("${pfm_dir}/*.pfm").collect()
        // Pair each .pfm file with the base name
        fils.collect { fil -> tuple(base, fil, kmerLen) }
        }.flatten().collate(3)
    
    // filter pfm_tuples to read-based and peak-based motifs
    pfm_tuples_reads = pfm_tuples.filter { base, pfm, kmerLen ->
        !base.contains("PEAKS")
    }
    pfm_tuples_peaks = pfm_tuples.filter { base, pfm, kmerLen ->
        base =~ /PEAKS/
    }

    // Combine motif results with prep output (pairing based on order of emission)
    pfm_bmark_prep_tuples_reads = pfm_tuples_reads.combine(prep, by: 0)
    pfm_bmark_prep_tuples_peaks = pfm_tuples_peaks.combine(prep_peaks, by: 0)
    // combine pfm_bmark_prep_top_frac_tuples_reads with a channel of top fractions (already done for peaks)
    pfm_bmark_prep_top_frac_tuples_reads = pfm_bmark_prep_tuples_reads.combine(top_fracs)
    // Run the benchmarking process with the combined motifs
    MOTIF_BENCHMARK(pfm_bmark_prep_top_frac_tuples_reads)
    MOTIF_BENCHMARK_PEAKS(pfm_bmark_prep_tuples_peaks)
    LOGO_MAKER(pfm_tuples)
    }