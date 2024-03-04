
# libraries
library(SingleCellExperiment)
library(tidyverse)
library(parallel)
library(scuttle)
library(MAST)

### data

sce_list <- list.files(
  pattern = ".*_sce_qced.rds",
  recursive = TRUE,
  full.names = TRUE
)
names(sce_list) <- c("LUAD", "K562")

sce_list <- lapply(sce_list, function(x){
  y <- readRDS(x)
  y[,y$sample %in% c(
    "K562-09", "K562-12",
    "LUAD-02", "LUAD-03"
  )]
})

metadata <- tibble(
  sample = c(
    "K562-09", "K562-12", "LUAD-02", "LUAD-03"
  ),
  type = rep(c("IVT", "TS"),2)
)

# add meta info
sce_list <- lapply(sce_list, function(x){
  colData(x) <- merge(colData(x), metadata)
  return(x)
})

### diff gene expression using MAST

## scale gene detection rate
sce_list <- mclapply(sce_list, function(x){
  colData(x)$n_genes <- scale(colData(x)$detected)
  x <- SceToSingleCellAssay(x, check_sanity = FALSE)
  return(x)
}, mc.cores	= 4, mc.preschedule = FALSE)

# assign reference group
sce_list <- lapply(sce_list, function(x){
  colData(x)$type <- as.factor(colData(x)$type)
  colData(x)$type <- relevel(colData(x)$type, "IVT")
  return(x)
})

# filter out non-expressed genes
sce_list <- lapply(sce_list, function(x){
  x <- x[rowSums(assay(x)) != 0, ]
  return(x)
})

# define & run hurdle model
zlm_list <- lapply(sce_list, function(x){
  zlm(formula = ~type + n_genes, sca = x, exprs_value = 'sf_lognorm')
  })

# get estimates
summary_zlm_list <- lapply(zlm_list, function(x){
    summary(x, doLRT= paste0("type", "TS"))
})

# get results
summary_zlm_tb_list <- lapply(summary_zlm_list, function(x){
  merge(
    x$datatable[contrast == "typeTS" & component=="H",.(primerid, `Pr(>Chisq)`)],
    x$datatable[contrast == "typeTS" & component=="logFC", .(primerid, coef)],
    by= "primerid") %>%
    mutate(
      FDR = p.adjust(`Pr(>Chisq)`, 'fdr'),
      primerid = gsub("GRCh38_", "", primerid)
    ) %>% 
    rename(
      gene_id = "primerid"
    )
})

# save
saveRDS(
  summary_zlm_tb_list,
  file="summary_zlm_tb_list.rds"
  )
