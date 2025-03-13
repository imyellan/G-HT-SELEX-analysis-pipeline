
process RUN_BOWTIE {
    conda "${HOME}/micromamba"
    label 'multithread'
    input:
        tuple val(baseName), path(pair1), path(pair2)

    output:
        tuple path("mapped_reads.bam"), val(baseName)

    script:
        """
        bowtie2 -p 4 --very-sensitive --no-unal --no-discordant --no-mixed \
        -x /home/hughespub/SELEX_Data/RawFlowcell_files/BowtieIndex/hg38 \
        -1 $pair1 -2 $pair2 > mapped_reads.sam
        samtools view -b -h mapped_reads.sam | samtools sort -o mapped_reads.bam
        rm mapped_reads.sam
        """
}

process JOIN_GHT_READS {
    conda "${HOME}/micromamba"
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

process PEAK_CALL {
    conda "${HOME}/micromamba"
    label 'high_mem'
    input:
        tuple val(baseName), path(bams), val(controlBase), path(controlBams)

    output:
        tuple path("macs_peaks.narrowPeak"), val(baseName), val("macs")
    
    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        bams_full = ${bams.join(' ')}
        controls_full = ${controlBams.join(' ')}
        """
        macs3 callpeak \
        -t $bams_full \
        -c $controls_full \
        -f BAMPE \
        -g hs \
        -n macs \
        -B -q 0.01
        """
}

process KNEEDLE_PEAKS {
    conda '/home/hugheslab1/iyellan/micromamba'
    input:
        tuple path(bed), val(baseName), val(peakSetType)

    output:
        tuple path("knee_filt_peaks.bed"), val(new_base), val("kneedle")

    publishDir "${params.outputDir}/${new_base}", mode: 'copy', overwrite: true

    script:
        new_base = baseName + "_kneedle"
        """
        magix_knee_filter.py $bed
        """
}

process SPLIT_PEAKS {
    conda "${HOME}/micromamba"
    input:
        tuple path(bed), val(baseName)

    output:
        tuple path("train.bed"), val(baseName), emit: train_bed
        tuple path("test.bed"), val(baseName), emit: test_bed
    
    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true
    publishDir (
        path: "${params.outputDir}/${baseName}",
        mode: 'copy',
        overwrite: true
    )

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
    conda "${HOME}/micromamba"
    //nPeaks in the output can serve as a placeholder instead of kmerLen
    input:
        tuple path(shuff_bed), val(baseName), val(nPeaks)

    output:
        tuple path("top_${nPeaks}.fasta"), val(baseName), val(nPeaks)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true, saveAs: {fn ->
        "${nPeaks}" + "_" + "$fn"
    }

    script:
        """
        if [[ $nPeaks -gt \$(cat $shuff_bed | wc -l) ]]; then
            # if nPeaks is greater than the number of peaks in the bed file, use the total # of peaks in the bed file
            num_peaks=\$(cat $shuff_bed | wc -l)
        else
            num_peaks=$nPeaks
        fi
        sort -k5 -gr $shuff_bed | head -n \$num_peaks \
        | awk 'BEGIN{OFS=FS="\t"} {print \$1,\$2,\$3,\$4,\$5}' > top_${nPeaks}.bed
        bedtools getfasta -fi \$HOME/data/hg38.fa \
        -bed top_${nPeaks}.bed -fo top_${nPeaks}.fasta
        """
}