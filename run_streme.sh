#!/bin/bash

DIR=$1
# Loop through all files ending with .fa in the current directory
for file in "$DIR"/*.fa; do
    # Extract the base name without the extension
    base_name="${file%.fa}"
    
    # Define the output directory based on the base name
    output_dir="${base_name}_streme"
    
    # Submit the job using submitjob with the specified options
    submitjob -w 60 -m 45 streme --p "$file" --thresh 0.0005 --oc "$output_dir" --minw 5 --maxw 15 --neval 40 --nref 8 --niter 25
done


# 477466.rc01
# 477467.rc01
# 477468.rc01
# 477469.rc01
# 477470.rc01
# 477471.rc01
# 477472.rc01
# 477473.rc01
# 477474.rc01
# 477475.rc01
# 477476.rc01
# 477477.rc01
# 477478.rc01
# 477479.rc01
# 477480.rc01
# 477481.rc01
# 477482.rc01
# 477483.rc01 meme homer dimant 