#!/bin/bash

DIR=$1

NUM_PROC=$(nproc)
# Loop through all files ending with .fa in the current directory
for file in "$DIR"/*.fa; do
    # Extract the base name without the extension
    base_name="${file%.fa}"
    
    # Define the output directory based on the base name
    output_dir="${base_name}_homer"
    # mem should be 15G if 20k in filename, else 30G
    # w should be 36H if 20k in filename, else 72H
    if [[ $file == *"20k"* ]]; then
        mem="15"
        w="36"
    else
        mem="30"
        w="72"
    fi
    # fasta-shuffle-letters -kmer 2 $file to temp file 
    dummy_file="/home/hugheslab1/spour/isaac_2024/tmp/$(basename "$file")_shuffled.fa"
    fasta-shuffle-letters -kmer 2 "$file" > $dummy_file

    
    # Submit the job using submitjob with the specified options

    submitjob -w $w -m $mem homer2 denovo -i "$file"  -p $NUM_PROC -b $dummy_file \> $output_dir
done

# 477628.rc01
# 477629.rc01
# 477630.rc01
# 477631.rc01
# 477632.rc01
# 477633.rc01
# 477634.rc01
# 477635.rc01
# 477636.rc01
# 477637.rc01
# 477638.rc01
# 477639.rc01
# 477640.rc01
# 477641.rc01
# 477642.rc01
# 477643.rc01
# 477644.rc01
# 477645.rc01