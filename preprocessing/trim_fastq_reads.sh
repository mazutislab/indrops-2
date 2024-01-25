#!/bin/bash

# activate conda-env
conda  activate indrops-2

# get files to array
fastqs=( $(find $PWD -type f -name "*R2_001.fastq.gz") )

#run loop
for fastq in ${fastqs[@]}; do
    
    # get directory and basename
    prefix="${fastq%%.*}"
    basename=$(basename "${fastq}")
    
    echo 'processing: ' ${fastq}
        
    zcat ${fastq} \
    | seqkit subseq --region 1:35 \
    -o ${prefix}/trimmed/${basename}
    
    echo 'writting to: ' ${prefix}/trimmed/${basename}

done

# deactivate env
conda deactivate
