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
    mkdir -p ${prefix}/dowsampling_by_frac
    
    echo 'processing: ' ${bam}
    
    python scrna-utils/src/dowsample_cb_bam.py \
    --bam ${bam} \
    --rd_out ${prefix}/dowsampling_by_frac/read_out.tsv \
    --ub_out ${prefix}/dowsampling_by_frac/umis_out.tsv \
    --gn_out ${prefix}/dowsampling_by_frac/gene_out.tsv \
    --n 18
        
done

# deactivate env
conda deactivate
