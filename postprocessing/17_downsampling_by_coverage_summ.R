# libraries
# libraries
library(tidyverse)
library(Matrix)
library(purrr)

# dir containing output of downsample_by_coverage.sh
path <- "/path/to/K562/dowsampling_by_cov"

## downsampled reads

# get files
read_files <- list.files(
  path = path,
  pattern = "reads_per_cb.tsv",
  recursive = TRUE,
  full.names = TRUE
)

# name files
names(read_files) <- gsub(
  "/d", "", str_extract(read_files, "K562_.*/d")
)

## downsampled cell x gene

# get files get files
mtx_files <- list.files(
  path = path,
  pattern = "ub_mat.mtx",
  recursive = TRUE,
  full.names = TRUE
)

# name files
names(mtx_files) <- gsub(
  "/d", "", str_extract(mtx_files, "K562_.*/d")
)

# gene names
feature_files <- list.files(
  path = path,
  pattern = "ub_features.tsv",
  recursive = TRUE,
  full.names = TRUE
)

# name files
names(feature_files) <- gsub(
  "/d", "", str_extract(feature_files, "K562_.*/d")
)

# barcodes
barcode_files <- list.files(
  path = path,
  pattern = "ub_barcodes.tsv",
  recursive = TRUE,
  full.names = TRUE
)

# name files
names(barcode_files) <- gsub(
  "/d", "", str_extract(barcode_files, "K562_.*/d")
)

# count matrix
count_list <- lapply(mtx_files, readMM)

# read
feature_list <- lapply(feature_files, function(x){
  x %>% read_tsv(col_names = F) %>% pull(X1)
})
barcode_list <- lapply(barcode_files, function(x){
  x %>% read_tsv(col_names = F) %>% pull(X1)
})

### combine
counts_list <- lapply(seq(count_list), function(i){
  colnames(count_list[[i]]) <- feature_list[[i]]
  rownames(count_list[[i]]) <- barcode_list[[i]]
  return(count_list[[i]])
})
names(counts_list) <- names(count_list)

# remove all zero
counts_list <- lapply(counts_list, function(x){
  idx <- rowSums(x) !=0
  idx <- idx[!is.na(idx)]
  m <- x[idx,]
  return(m)
})

# umi counts
counts_umi <- lapply(counts_list, function(x){
  y <- rowSums(x)
  y %>% stack() %>% as_tibble() %>% 
    dplyr::rename(
      "barcode" = "ind", "umis" = "values"
    )
  }) %>% 
  bind_rows(
    .id = "library"
  ) %>% 
  filter(!barcode %in% "-")

# gene counts
counts_gene <- lapply(counts_list, function(x){
  m <- x %>% as(., "nsparseMatrix")*1
  m %>% rowSums() %>% 
    as_tibble(rownames = "barcode") %>% 
    dplyr::rename(gene = "value")
  }) %>% 
  bind_rows(
    .id = "library"
  ) %>% 
  filter(!barcode %in% "-")

# read counts
read_counts <- lapply(read_files, read_tsv) %>% 
  bind_rows(
    .id = "library"
  )  %>% 
  filter(!barcode %in% "-")

# combine data
counts <- left_join(
  read_counts, counts_umi
  ) %>% 
  left_join(
    counts_gene
  )

write_tsv(
  counts, 
  file = "K562_counts_to_20k_reads.tsv.gz"
)
