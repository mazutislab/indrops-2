# libraries
library(tidySingleCellExperiment)
library(SingleCellExperiment)
library(tidyverse)
library(scuttle)

# read
sce_qced <- readRDS(
  file = "K562_sce_qced.rds"
)

# frac per cell
gene_len_cat_tb <- sce_qced %>%
  select(
    sample, matches("subsets_[0-9]{1,2}gtile_percent")
  ) %>% 
  gather(variable, value, -sample) %>% 
  mutate(
    variable = gsub("subsets_|gtile_percent", "", variable),
    variable = factor(variable, levels = 1:10)
  )

# write
write_tsv(
  gene_len_cat_tb,
  file = "K562_gene_length_cat_frac_per_cell.tsv.gz"
)

# avg expression per library
sce_qced_list <- lapply(unique(sce_qced$sample), function(x){
  sce_qced[, sce_qced$sample == x]
})
names(sce_qced_list) <- unique(sce_qced$sample)

len_lognorm <- lapply(sce_qced_list, function(x){
  assay(x, "sf_lognorm") %>% rowMeans2(useNames=TRUE) %>%
    stack() %>% as_tibble() %>% 
    rename("gene_id" = "ind") %>% 
    left_join(
      rowData(x) %>% as_tibble()
    )
  }) %>% 
  bind_rows(.id = "sample") %>% 
  mutate(
    gene_length_cat = as.factor(ntile(gene_length, n = 3))
  )

# write
write_tsv(
  len_cpm,
  file = "K562_len_lognorm_per_library.tsv.gz"
)
