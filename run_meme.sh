#!/bin/bash

DIR=$1

NUM_PROC=$(nproc)
# Loop through all files ending with .fa in the current directory
for file in "$DIR"/*.fa; do
    # Extract the base name without the extension
    base_name="${file%.fa}"
    
    # Define the output directory based on the base name
    output_dir="${base_name}_meme"
    
        # mem should be 15G if 20k in filename, else 30G
    # w should be 36H if 20k in filename, else 72H
    if [[ $file == *"20k"* ]]; then
        mem="15"
        w="36"
    else
        mem="40"
        w="72"
    fi
    
    # Submit the job using submitjob with the specified options
    submitjob -w $w -m $mem meme "$file"  -oc "$output_dir" -minw 5 -maxw 15 -p $NUM_PROC -hsfrac 0.1

done

# 477567.rc01
# 477568.rc01
# 477569.rc01
# 477570.rc01
# 477571.rc01
# 477572.rc01
# 477573.rc01
# 477574.rc01
# 477575.rc01
# 477576.rc01
# 477577.rc01
# 477578.rc01
# 477579.rc01
# 477580.rc01
# 477581.rc01
# 477582.rc01
# 477583.rc01
# 477584.rc01