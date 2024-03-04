# libraries
library(tidyverse)

# dir containing output of downsample_by_fraction.sh
path <- "/path/to/K562/dowsampling_by_frac"

# get files
read_list <- list.files(
  path = path,
  pattern = "read_out.tsv",
  recursive = TRUE,
  full.names = TRUE
)
names(read_list) <- gsub(
  "/d", "", str_extract(read_list, "K562_.*/d")
)

umis_list <- list.files(
  path = path,
  pattern = "umis_out.tsv",
  recursive = TRUE,
  full.names = TRUE
)
names(umis_list) <- gsub(
  "/d", "", str_extract(umis_list, "K562_.*/d")
)

gene_list <- list.files(
  path = path,
  pattern = "gene_out.tsv",
  recursive = TRUE,
  full.names = TRUE
)
names(gene_list) <- gsub(
  "/d", "", str_extract(gene_list, "K562_.*/d")
)

# reads 
read_tb <- lapply(read_list, read_tsv) %>% 
  bind_rows(.id = "sample") %>%
  gather(prop, reads, -sample, -barcode)

# umis
umis_tb <- lapply(umis_list, read_tsv) %>% 
  bind_rows(.id = "sample") %>%
  filter(`1.0` >= 1000) %>% 
  filter(!barcode %in% "-") %>% 
  group_by(sample) %>% 
  slice_max(order_by = `1.0`, prop=0.25) %>% 
  gather(prop, umis, -sample, -barcode) %>% 
  ungroup()

# genes
gene_tb <- lapply(gene_list, read_tsv) %>% 
  bind_rows(.id = "sample") %>%
  gather(prop, genes, -sample, -barcode)

# combine data
combined_tb <- read_tb %>% 
  right_join(., umis_tb) %>% 
  left_join(., gene_tb)

# summarize stats
summ_tb <- combined_tb %>%
  group_by(sample, prop) %>% 
  summarise_if(
    is.numeric,
    funs(mean, median, sum)
  ) %>% 
  ungroup()

# write output
write_tsv(
  summ_tb, 
  file = "K562_saturation.tsv.gz"
)

