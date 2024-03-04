#!/bin/bash

# activate conda-env
conda  activate indrops-2

# get reference 
human_bed_url="https://altushost-swe.dl.sourceforge.net/project/rseqc/BED/Human_Homo_sapiens/hg38_GENCODE_V42_Comprehensive.bed.gz"
human_bed="hg38_GENCODE_V42_Comprehensive.bed"

if [ ! -f "$human_bed" ]; then
    curl -sS "$human_bed_url" | zcat > "$human_bed"
fi

# get files to array
bams=( $(find $PWD/starsolo -type f -name "*Aligned.sortedByCoord.out.bam") )

#run loop
for bam in ${bams[@]}; do
    
    # get parent directory
    prefix=$(dirname $(dirname "${bam}"))
    
    # make directory for output
    mkdir -p ${prefix}/rseqc_stats
    
    echo 'processing: ' ${bam}
    
    geneBody_coverage.py \
        --input ${bam} \
        --refgene $human_bed \
        --out-prefix ${prefix}/rseqc_stats/
        
    read_distribution.py \
        --input ${bam} \
        --refgene $human_bed \
        > ${prefix}/rseqc_stats/read_distribution.txt
    
    read_GC.py \
        --input ${bam} \
        --out-prefix ${prefix}/rseqc_stats/

done

# deactivate env
conda deactivate