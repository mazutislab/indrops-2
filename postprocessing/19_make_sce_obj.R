
# libraries
library(tidySingleCellExperiment)
library(SingleCellExperiment)
library(tidyverse)
library(scuttle)

# path to data
path <- "/path/to/K562/starsolo/"

# get matrices
file_list <- list.files(
  path = path,
  pattern = "matrix.mtx|^features.tsv|^barcodes.tsv",
  recursive = TRUE,
  full.names = TRUE
)

# use solo filtered counts
file_list_filt <- file_list[
  grepl("filtered", file_list)
]

# name files
names(file_list_filt) <- gsub(
  "/s", "", str_extract(file_list_filt, "K562_.*/s")
)

# read gene info
gene_info <- read_tsv(
  "GRCh38_ensembl_v107_gencode_v41_gene_info.tsv"
)

# read to sce obj
sce_filt_list <- lapply(unique(names(file_list_filt)), function(s){
  # read mtx
  mat <- grep(paste0(s, "/", ".*matrix.mtx"), file_list_filt, value = TRUE)
  mat <- Matrix::readMM(mat)
  # read barcodes
  col <- grep(paste0(s, "/", ".*barcodes.tsv"), file_list_filt, value = TRUE)
  col <- read_tsv(col, col_names = FALSE)
  # read features
  row <- grep(paste0(s, "/", ".*features.tsv"), file_list_filt, value = TRUE)
  row <- read_tsv(row, col_names = FALSE) %>% 
    dplyr::rename("gene_id"="X1") %>% 
    select(gene_id)
  # assign names
  rownames(mat) <- row$gene_id
  colnames(mat) <- col$X1
  # make sce obj
  rowData <- DataFrame(
    left_join(row, gene_info)
  )
  sce <- SingleCellExperiment(
    assays = list(counts = mat),
    rowData = rowData
  )
  sce$sample <- s
  return(sce)
})
names(sce_filt_list) <- unique(names(file_list_filt))

# concatenate sce obj
sce_filt <- Reduce(cbind, sce_filt_list)

# remove zeros
sce_qced <- sce[which(rowSums2(counts(sce)) > 0)]

# keep cells with at least 1k UMIs
sce_qced <- sce_qced[,which(colSums2(counts(sce_qced)) >= 1000)]

# mitochondrial genes
mt <- rownames(sce_qced)[grep("^MT-", rowData(sce_qced)$gene_name)]

# ribosomal genes
rb <- rownames(sce_qced)[grep("^RP[SL]", rowData(sce_qced)$gene_name)]

# make a list of gene length categories
ll <- lapply(unique(rowData(sce_qced)$gene_length_cat), function(x){
  rownames(sce_qced)[grep(x, rowData(sce_qced)$gene_length_cat)]
})
names(ll) <- unique(rowData(sce_qced)$gene_length_cat)

# make a list of gene biotype categories
bl <- lapply(unique(rowData(sce_qced)$gene_type), function(x){
  rownames(sce_qced)[grep(x, rowData(sce_qced)$gene_type)]
})
names(bl) <- unique(rowData(sce_qced)$gene_type)

# calculate cq metrics
sce_qced <- addPerCellQC(
  sce_qced, flatten = T,
  subsets = c(
    list(mt = mt, ribo = rb), ll, bl
    )
)

# log normalize
assay(sce_qced, "sf_lognorm") <- normalizeCounts(sce_qced)

# cpm normalize
assay(sce_qced, "cpm") <- calculateCPM(sce_qced)

# save
saveRDS(
  sce_qced,
  file = "K562_sce_qced.rds"
)


