
# libraries
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)


# read annotations
gtf_file <- paste0(
  "gencode.v41.primary_assembly.annotation.gtf.filtered"
)

gtf <- import(
  gtf_file, format = "gtf"
)

# get genes
gene_info <- gtf %>% 
  as_tibble() %>%
  select(
    seqnames:width, type, gene_id, gene_type:gene_name
  ) %>% 
  filter(
    type %in% "gene"
  ) %>% 
  mutate(
    gene_length = width,
    gene_name = gsub("GRCh38_", "", gene_name),
    chr = gsub("GRCh38_", "", seqnames)
  ) %>% 
  select(
    gene_id:chr
  ) %>% 
  mutate(
    gene_length_cat = paste0(ntile(gene_length, 10),"tile")
  )

# write
write_tsv(
  gene_info, "GRCh38_ensembl_v107_gencode_v41_gene_info.tsv"
)