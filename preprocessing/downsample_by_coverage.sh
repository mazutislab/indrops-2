#!/bin/bash

# activate conda-env
conda  activate indrops-2

# get scripts
git clone https://github.com/jsimonas/scrna-utils

# get files to array
bams=( $(find $PWD/starsolo -type f -name "*Aligned.sortedByCoord.out.bam") )

#run loop
for bam in ${bams[@]}; do
    
    # get parent directory
    prefix=$(dirname $(dirname "${bam}"))
    
    # make directory for output
    mkdir -p ${prefix}/dowsampling_by_cov
    
    echo 'processing: ' ${bam}
    
    python scrna-utils/src/dowsample_bam_by_cb_coverage.py \
    --bam ${bam} \
    --out ${prefix}/dowsampling_by_cov \
    --target_num_reads 20000 \
    --n 18
        
done

# deactivate env
conda deactivate