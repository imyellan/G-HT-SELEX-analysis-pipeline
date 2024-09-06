#!/usr/bin/env nextflow

params.inputDir = "${HOME}/nextflow_test"
params.outputDir = "${HOME}/TEHMM_proj/analysis/selex_results/nextflow_test"

// Channel to read all fastq files in the input directory
Channel
    .fromPath("${params.inputDir}/*.fastq.gz")
    .map {fastq ->
        def baseName = fastq.baseName.replaceAll(/\.fastq\.gz$/, '')
        tuple(fastq, baseName)
    }
    .set { fastq_tuples }

// print the project directory
println "Project directory: ${projectDir}"

// Define the process to deduplicate the fastq files
process DE_DUP {
    conda '/home/hugheslab1/iyellan/micromamba'
    input:
        tuple path(fastq), val(baseName)

    output:
        tuple path("dedupe.fastq.gz"), val(baseName)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy'

    script:
        """
        dedupe.sh in=$fastq out="dedupe.fastq.gz"
        """
}

// Use reformat.sh to convert fastq to fasta and split into 70% train and 30% test
process SPLIT {
    conda '/home/hugheslab1/iyellan/micromamba'
    input:
        tuple path(fastq), val(baseName)

    output:
        tuple path("train.fastq.gz"), val(baseName), emit: train
        tuple path("test.fastq.gz"), val(baseName), emit: test

    publishDir "${params.outputDir}/${baseName}", mode: 'copy'

    script:
        """
        seqkit shuffle --rand-seed 1 $fastq -o shuffled.fastq.gz
        seqkit sample -p 0.3 shuffled.fastq.gz -o test.fastq.gz
        seqkit seq -n -i test.fastq.gz -o test_ids.txt
        seqkit grep -v -f test_ids.txt shuffled.fastq.gz -o train.fastq.gz
        """
}

// Use the autosome-ru folk's k-mer enrichment script to subset the training data fastq 
// to the 10,000 sequences with the highest k-mer enrichment scores

process KMER_ENRICH {
    conda '/home/hugheslab1/iyellan/micromamba'
    input:
        tuple path(fastq), val(baseName)

    output:
        tuple path("kmer_enr.fasta.gz"), val(baseName), emit: kmerEnr

    publishDir "${params.outputDir}/${baseName}", mode: 'copy'

    script:
        """
        $HOME/software/HT-SELEX-kmer-filtering/extract_topk.py -d -k 5 \
        --fastq $fastq --out kmer_enr.fasta.gz -b 10000 -a 0
        """
}

process SEQ_PREP {
    conda '/home/hugheslab1/iyellan/micromamba'
    input:
        tuple path(fasta), val(baseName)

    output:
        // tuple path("unzipped.fasta"), val(baseName), emit: unzipped
        tuple path("unzipped_unique_ids.fasta"), val(baseName), emit: unzipped

    // publishDir "${params.outputDir}/${baseName}", mode: 'copy'

    script:
        """
        gunzip "$fasta" -c > unzipped.fasta
        cat unzipped.fasta | fasta-unique-names > unzipped_unique_ids.fasta
        """

}

process MOTIF_DISCOVERY {
    label 'multithread'
    conda '/home/hugheslab1/iyellan/micromamba'

    input:
        tuple path(fasta), val(baseName)

    output:
        tuple path("meme_out/*"), val(baseName), emit: meme
        tuple path("streme_out/*"), val(baseName), emit: streme
        tuple path("homer_out.txt"), val(baseName), emit: homer

    publishDir "${params.outputDir}/${baseName}", mode: 'copy'

    script:
        """
        ## MEME
        meme $fasta -oc meme_out/ -minw 5 -maxw 15 -p 4
        ## STREME
        streme --p $fasta --oc streme_out/ --minw 5 --maxw 15 --neval 40 --nref 8 --niter 25
        ## HOMER
        fasta-shuffle-letters -kmer 2 -seed 42 $fasta > dummy_file
        homer2 denovo -i $fasta -p 4 -b dummy_file -o homer_out.txt
        ## DIMONT
        #java -jar dimont/Dimont.jar data=unzipped_unique.fasta infix=dimont_out 
        """
}

process RUN_BEESEM {
    conda '/home/hugheslab1/iyellan/micromamba/envs/py27_env'
    label 'RUN_BEESEM'
    input:
        tuple path(fasta), val(baseName)
    
    output:
        tuple path("beesem_out/"), val(baseName), emit: beesem
    
    publishDir "${params.outputDir}/${baseName}", mode: 'copy'

    script:
        """
        $HOME/scripts/fasta_linearize.sh $fasta \
        | awk -F"\t" '{print \$2}' | uniq -c | awk 'BEGIN{OFS="\t"} {print \$2,\$1}' \
        > beesem_formatted.fasta
        python $HOME/software/BEESEM/beesem.py -o beesem_out "${baseName}" beesem_formatted.fasta
        """
}

// run DE_DUP and SPLIT as workflow
workflow {
    deduped = DE_DUP(fastq_tuples)
    split = SPLIT(deduped)
    kmerEnr = KMER_ENRICH(split.train)
    unzipped = SEQ_PREP(kmerEnr)
    motifs = MOTIF_DISCOVERY(unzipped)
    beesem = RUN_BEESEM(unzipped)
}