process MOTIF_DISCOVERY {
    label 'multithread'
    conda "${HOME}/micromamba"

    input:
        tuple path(fasta), val(baseName), val(kmerLen), val(motifCaller)

    output:
        tuple path("${motifCaller}_out_${kmerLen}/"), val(baseName), val(kmerLen)

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        if("${motifCaller}" == "meme")
            """
            meme $fasta -oc meme_out_${kmerLen}/ -minw 5 -maxw 20 -p 4 -dna -nmotifs 3
            meme_xml_parse.py meme_out_${kmerLen}/meme.xml
            """
        
        else if("${motifCaller}" == "streme")
            """
            streme --p $fasta --oc streme_out_${kmerLen}/ --minw 5 --maxw 20 --neval 40 --nref 8 --niter 25 --dna
            meme_xml_parse.py streme_out_${kmerLen}/streme.xml
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
    conda "${HOME}/micromamba/envs/py27_env"
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

process RUN_DIMONT_GHT{
    label 'multithread_dimont'
    conda "${HOME}/micromamba"

    input:
        tuple path(peaks), val(baseName)

    output:
        tuple path("dimont_ght_out_0/"), val(baseName), env('NPeaks'), emit: dimont_ght

    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        libdir=/home/hugheslab1/iyellan/software/Jstacs/lib/xml-commons

        ## extract fastas in right format - also remove peaks with Ns
        awk 'BEGIN{OFS=FS="\t"} {print \$1,\$2,\$3," peak: 100; signal: "\$5}' $peaks \
        | bedtools getfasta -fi $HOME/data/hg38.fa -bed - -nameOnly \
        | ~/scripts/fasta_linearize.sh \
        | awk -F"\t" 'BEGIN{OFS="\\n"} \$2!~/[nN]/ {print "> "gensub(","," ","g",\$1),\$2}' \
        > peaks.fa
        NPeaks=\$(grep ">" peaks.fa | wc -l | awk '{print \$1}')

        java -Djava.awt.headless=true -cp \${libdir}/batik-anim.jar:\${libdir}/batik-awt-util.jar:\${libdir}/batik-bridge.jar:\${libdir}/batik-codec.jar:\${libdir}/batik-css.jar:\${libdir}/batik-dom.jar:\${libdir}/batik-extension.jar:\${libdir}/batik-ext.jar:\${libdir}/batik-gui-util.jar:\${libdir}/batik-gvt.jar:\${libdir}/batik.jar:\${libdir}/batik-parser.jar:\${libdir}/batik-script.jar:\${libdir}/batik-svg-dom.jar:\${libdir}/batik-svggen.jar:\${libdir}/batik-swing.jar:\${libdir}/batik-transcoder.jar:\${libdir}/batik-util.jar:\${libdir}/batik-xml.jar:\${libdir}/xml-apis-ext.jar:\${libdir}/pdf-transcoder.jar:$HOME/build/Dimont-HTS/Dimont-HTS.jar \
        projects.dimont.Dimont data=peaks.fa threads=4 infix=dimont_ght position=peak value=signal

        ## parse xml outputs to pfm using Ali's scripts and then my pwm to pfm script; 
        ## move results to single directory
        mkdir -p dimont_ght_out_0
        for xml_model in *_ght-pwm-*.xml; do
            model_num=\$(basename \${xml_model} | sed -E 's/.*([0-9]+).xml/\1/g')
            dimont_motif_parse.R \${xml_model}
            #DimontMotifParser.py \
            #-i \${xml_model} -o dimont_ght_out_0 -n dimont_ght_pfm_\${model_num}
            #reformat_dimont_ght_motif.R dimont_ght_out_0/dimont_ght_pfm_\${model_num}.pwm
        done
        rm peaks.fa
        mv *.xml *.png *.pfm dimont_ght_out_0
        """
}

process RUN_DIMONT_HTS{
    label 'multithread_dimont'
    conda "${HOME}/micromamba"

    input:
        tuple path(fasta), val(baseName)

    output:
        tuple path("dimont_hts_out_0/"), val(baseName), val(0), emit: dimont_hts
    
    publishDir "${params.outputDir}/${baseName}", mode: 'copy', overwrite: true

    script:
        """
        gunzip -c $fasta > in_fasta.fa
        libdir=/home/hugheslab1/iyellan/software/Jstacs/lib/xml-commons
        java -Djava.awt.headless=true -cp \${libdir}/batik-anim.jar:\${libdir}/batik-awt-util.jar:\${libdir}/batik-bridge.jar:\${libdir}/batik-codec.jar:\${libdir}/batik-css.jar:\${libdir}/batik-dom.jar:\${libdir}/batik-extension.jar:\${libdir}/batik-ext.jar:\${libdir}/batik-gui-util.jar:\${libdir}/batik-gvt.jar:\${libdir}/batik.jar:\${libdir}/batik-parser.jar:\${libdir}/batik-script.jar:\${libdir}/batik-svg-dom.jar:\${libdir}/batik-svggen.jar:\${libdir}/batik-swing.jar:\${libdir}/batik-transcoder.jar:\${libdir}/batik-util.jar:\${libdir}/batik-xml.jar:\${libdir}/xml-apis-ext.jar:\${libdir}/pdf-transcoder.jar:$HOME/build/Dimont-HTS/Dimont-HTS.jar \
        projects.dimont.hts.DimontTool i=in_fasta.fa n=4

        ## move results to single directory
        mkdir dimont_hts_out_0
        for motifdir in Motif_*; do
            motif_num=\$(basename \$motifdir | sed 's/Motif_//g')
            tail -n +2 \${motifdir}/Model_PWM.pwm | head -n -1 \
            > dimont_hts_out_0/dimont_pfm_\${motif_num}.pfm
            if [[ -f \${motifdir}/Sequence_logo_\${motif_num}.pdf ]]; then
                mv \${motifdir}/*xml \${motifdir}/*tsv \${motifdir}/*pdf dimont_hts_out_0
            fi
        done
        rm in_fasta.fa
        """
}