#3 Annotation

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse'))

# Annotation---------------------------------------------------------------
se <- readRDS(paste0(save.path,'postClust.rds'))

new.ident <- c(
  '0' = "ZF",
  '1' = "ZR",
  '2' = "Capsule", 
  '3' = "ZF",
  '4' = "ZG",
  '5' = 'ZR',
  '6' = 'ZF',
  '7' = 'ZF',
  '8' = 'Stromal',
  '9' = 'Medulla',
  '10' = 'ZG',
  '11' = 'Stromal'
)

Idents(se) <- se$seurat_clusters
se <- RenameIdents(se, new.ident)
se$annotation <- Idents(se) <- factor(Idents(se),levels = c("ZG","ZF","ZR","Capsule","Medulla",'Stromal'))

# save
saveRDS(se, paste0(save.path,'postAnnot.rds'))

# DE analysis---------------------------------------------------------------
Idents(se) <- se$annotation
DEG <- FindAllMarkers(se,assay = 'SCT',logfc.threshold = 0.25,densify = T)