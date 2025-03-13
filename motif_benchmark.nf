
process BENCHMARK_PREP{
    label = 'benchmark_reads'
    input:
    // tuple path(motif_dir), val(baseName)
        tuple path(test_fq), val(baseName)

    output:
        tuple val(baseName), path("pos.fa.gz"), path("neg.fa.gz")

    publishDir (
        path: "${params.outputDir}/${baseName}", 
        mode: 'copy', 
        overwrite: true
    )
    
    script:
    """
    flank_5=ACACTCTTTCCCTACACGACGCTCTTCCGATCT\$(echo "${baseName}" | sed -E 's/.*_([ATCG]{2})40N.*/\\1/g')
    flank_3=\$(echo "${baseName}" | sed -E 's/.*40N([ATCG]+)_.*/\\1/g')AGATCGGAAGAGCACACGTCTGAACTCCAG
    #prepare --seq $test_fq --positive-file pos.fa.gz --negative-file neg.fa.gz \
    #--flank-5 \$flank_5 --flank-3 \$flank_3 --seed 42
    prepare --seq $test_fq --positive-file pos.fa.gz --negative-file neg.fa.gz \
    --seed 42
    """
}

process BENCHMARK_PREP_PEAKS{
    label = 'benchmark_peaks'
    input:
    // tuple path(motif_dir), val(baseName)
        tuple path(test_bed), val(baseName), val(top_frac)

    output:
        tuple val(baseName), path("pos.fa.gz"), path("neg.fa.gz"), val(top_frac)

    publishDir (
        path: "${params.outputDir}/${baseName}", 
        mode: 'copy', 
        overwrite: true
    )
    
    script:
    """
    n_peaks=\$(cat ${test_bed} | wc -l)
    top_n_as_pos=\$(printf "%.0f" \$(echo "${top_frac}*\${n_peaks}" | bc))
    if [[ "\$top_n_as_pos" -eq 0 ]]; then
        top_n_as_pos=1
    fi
    awk 'BEGIN{OFS=FS="\t"} {print \$1, \$2, \$3, \$4, \$5}' ${test_bed} > test_tmp.bed
    prepare --peaks test_tmp.bed --bed --positive-file pos.fa.gz --negative-file neg.fa.gz \
    --seed 42 --assembly-fasta \$HOME/data/hg38.fa --assembly-sizes \$HOME/data/hg38.chrom.sizes \
    --top \$top_n_as_pos
    """
}

process MOTIF_BENCHMARK{
    label = 'benchmark_reads'
    time = '3:00:00'
    input:
        tuple val(baseName), path(motif), val(kmer_len), path(pos_fa), path(neg_fa), val(top_frac)

    // output motif_eval directory
    output:
        tuple path("ROC.tsv"), path("PR.tsv"), path("eval_summ.txt")

    publishDir (
        path: "${params.outputDir}/${baseName}",
        mode: 'copy',
        overwrite: true,
        saveAs: {fn ->
            def motifbase = motif.getBaseName()
            // remove periods from top_frac
            def top_frac_no_dot = top_frac.toString().replaceAll(/\./, "")
            "$motifbase" + "_" + "$kmer_len" + "_" + "$top_frac_no_dot" + "_" + "$fn"
        }
    )
    
    script:
    """
    evaluate --positive-file $pos_fa --negative-file $neg_fa --top $top_frac \
    --motif $motif --roc --pr --pfm --roc-filename ROC.tsv --pr-filename PR.tsv \
    --bins 10000 > eval_summ.txt
    #--plot --plot-filename ROC.png --json
    """
}

process MOTIF_BENCHMARK_PEAKS{
    label = 'benchmark_peaks'
    time = '3:00:00'
    input:
        tuple val(baseName), path(motif), val(kmer_len), path(pos_fa), path(neg_fa), val(top_frac)

    // output motif_eval directory
    output:
        tuple path("ROC.tsv"), path("eval_summ.txt")

    publishDir (
        path: "${params.outputDir}/${baseName}",
        mode: 'copy',
        overwrite: true,
        saveAs: {fn ->
            def motifbase = motif.getBaseName()
            // remove periods from top_frac
            def top_frac_no_dot = top_frac.toString().replaceAll(/\./, "")
            "$motifbase" + "_" + "$kmer_len" + "_" + "$top_frac_no_dot" + "_" + "$fn"
        }
    )
    
    script:
    """
    evaluate --positive-file $pos_fa --negative-file $neg_fa \
    --motif $motif --roc --pfm --roc-filename ROC.tsv \
    --bins 10000 > eval_summ.txt
    """
}

process LOGO_MAKER{
    conda '/home/hugheslab1/iyellan/micromamba'
    label = 'logo'

    input:
        tuple val(baseName), path(motif), val(kmer_len)

    output:
        tuple path("logo.png"), path("logo_REV.png")

    publishDir (
        path: "${params.outputDir}/${baseName}",
        mode: 'copy',
        saveAs: {fn ->
            def motifbase = motif.getBaseName()
            "$motifbase" + "_" + "$kmer_len" + "_" + "$fn"
        },
        overwrite: true
    )

    script:
    """
        plot_motif_logo_kaitlin_mod.R $motif
    """
}


