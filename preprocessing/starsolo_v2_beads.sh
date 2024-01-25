#!/bin/bash

# activate conda-env
conda  activate indrops-2

# clone repo to get barcodes
git clone https://github.com/mazutislab/indrops-2

# barcode white list  
bc1=( $(find $PWD -type f -name "v2_cb1.txt") )
bc2=( $(find $PWD -type f -name "v2_cb2.txt") )

# get barcode dirs to array
barcode=( $(find $PWD -type f -name "*R1_001.fastq.gz") )

# get genomic dirs to array
genomic=( $(find $PWD -type f -name "*R2_001.fastq.gz") )

# run starsolo
for i in "${!genomic[@]}"; do
    
    # get file name
    prefix="${genomic[$i]%%.*}"
    
    echo 'processing: ' $(basename "${prefix}")
    echo 'genomic: ' ${genomic[$i]}
    echo 'barcode: ' ${barcode[$i]}

    STAR \
        --genomeDir 'path/to/genome/index' \
        --readFilesIn ${genomic[$i]} ${barcode[$i]} \
        --soloCBwhitelist $bc1 $bc2 \
        --runThreadN 16 \
        --outFileNamePrefix ${prefix}/starsolo/ \
        --sjdbOverhang 100 \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingBinsN 20 \
        --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
        --runDirPerm All_RWX \
        --readFilesCommand zcat \
        --soloFeatures GeneFull \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype EditDist_2 \
        --soloUMIdedup Exact \
        --soloCBposition 0_0_0_7 0_12_0_19 \
        --soloUMIposition 0_20_0_27
    
    echo 'done'

done


# deactivate env
conda deactivate
