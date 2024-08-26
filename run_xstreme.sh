#!/bin/bash

DIR=$1

NUM_PROC=$(nproc)
# Loop through all files ending with .fa in the current directory
for file in "$DIR"/*.fa; do
    # Extract the base name without the extension
    base_name="${file%.fa}"
    
    # Define the output directory based on the base name
    output_dir="${base_name}_xstreme"
    
    # Submit the job using submitjob with the specified options
    submitjob -w 60 -m 45 xstreme --p "$file"  --oc "$output_dir" --minw 5 --maxw 15 --meme-p $NUM_PROC 
done

# 477514.rc01
# 477515.rc01
# 477516.rc01
# 477517.rc01
# 477518.rc01
# 477519.rc01
# 477520.rc01
# 477521.rc01
# 477522.rc01
# 477523.rc01
# 477524.rc01
# 477525.rc01
# 477526.rc01
# 477527.rc01
# 477528.rc01
# 477529.rc01
# 477530.rc01
# 477531.rc01