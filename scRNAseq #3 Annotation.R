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
  "0" = "ZF_3",
  "1" = "ZF_2",
  "2" = "Macrophage",
  "3" = "Fibroblast",
  "4" = "Capsule",
  "5" = "ZF_1",
  "6" = "Endothelial",
  "7" = "ZR",
  "8" = "Medulla",
  "9" = "Tcell",
  "10" = "VSM",
  "11" = "ZG",
  "12" = "Adipocyte",
  "13" = "Bcell",
  "14" = "ZF_4",
  "15" = "Neutrophil"
)
se$sub_annotation <- new.ident[se$seurat_clusters]
se$sub_annotation <- factor(se$sub_annotation, 
                        levels = c("ZG","ZF_1","ZF_2","ZF_3","ZF_4","ZR",
                                   "Capsule","Medulla","Endothelial",
                                   'Fibroblast',"VSM",'Adipocyte',
                                   "Macrophage","Tcell","Bcell",'Neutrophil'))

se$annotation <- sapply(strsplit(new.ident,split = '_'),'[[',1)[se$seurat_clusters]
se$annotation <- factor(se$annotation, 
                        levels = c("ZG","ZF","ZR","Capsule","Medulla",
                                   "Endothelial",'Fibroblast',"VSM",'Adipocyte',
                                   "Macrophage","Tcell","Bcell",'Neutrophil'))

# save
saveRDS(se, paste0(save.path,'postAnnot.rds'))

# DE analysis---------------------------------------------------------------
Idents(se) <- se$annotation
DEG <- FindAllMarkers(se,assay = 'SCT',logfc.threshold = 0.25,densify = T)