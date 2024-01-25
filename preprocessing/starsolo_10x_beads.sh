#!/bin/bash

# activate conda-env
conda  activate indrops-2

# clone repo to get barcodes
git clone https://github.com/mazutislab/indrops-2

# barcode white list  
bc_whitelist='path/to/3M-february-2018.txt'

# get barcode dirs to array
barcode=( $(find $PWD -type f -name "*R1_001.fastq.gz") )

# get genomic dirs to array
genomic=( $(find $PWD -type f -name "*R2_001.fastq.gz") )

# run starsolo
for i in "${!genomic[@]}"; do
    
    # get parent directory
    prefix= "${genomic[$i]%%.*}"
    
    echo 'processing: ' $(basename "${prefix}")
    echo 'genomic: ' ${genomic[$i]}
    echo 'barcode: ' ${barcode[$i]}

    STAR \
        --genomeDir 'path/to/genome/index' \
        --readFilesIn ${genomic[$i]} ${barcode[$i]} \
        --soloCBwhitelist $bc_whitelist \
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
        --soloType CB_UMI_Simple \
        --soloBarcodeReadLength 0 \
        --soloCBmatchWLtype 1MM_multi \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR
    
    echo 'done'

done


# deactivate env
conda deactivate
