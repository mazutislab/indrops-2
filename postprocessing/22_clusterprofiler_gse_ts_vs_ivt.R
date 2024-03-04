
# libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

### data

# save
summary_zlm_tb_list <- readRDS(
  file="summary_zlm_tb_list.rds"
)

# get entrez IDs for universe
ensembl <- useEnsembl(
  biomart = "genes"
)

ensembl <- useDataset(
  dataset = "hsapiens_gene_ensembl",
  mart = ensembl
)

universe_list <- lapply(summary_zlm_tb_list, function(x){
  getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
    filters="ensembl_gene_id",
    values = x$gene_id,
    mart=ensembl
  ) %>% 
    dplyr::rename(gene_id = "ensembl_gene_id")
})

dge_lfc_list <- lapply(summary_zlm_tb_list, function(x){
  l <- x$coef
  names(l) <- x$gene_id
  l <- na.omit(l)
  l = sort(l, decreasing = TRUE)
  l
})

gse_go_list <- lapply(dge_lfc_list, function(x){
  gseGO(
    geneList = x,
    ont ="ALL",
    keyType = "ENSEMBL", 
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.05,
    eps = 0,
    nPermSimple = 100000,
    verbose = TRUE,
    OrgDb = org.Hs.eg.db,
    pAdjustMethod = "none"
  )
})

# flatten the table
gse_go_list_tb <-gse_go_list %>% 
  lapply(., function(x){
    x %>% as_tibble() %>% 
      filter(p.adjust<0.05) %>%
      mutate(
        sign = case_when(
          enrichmentScore > 0 ~ "enriched in TS",
          TRUE ~ "enriched in IVT"
        )
      )
  }) %>% bind_rows(
    .id = "source"
  )

# write
write_tsv(
  gse_go_list_tb,
  "gse_go_tb.tsv"
)


