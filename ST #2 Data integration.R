#2 Normalization, Dimention reduction, Batch correction

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','harmony'))

# SCTransform---------------------------------------------------------------
se.list <- readRDS(paste0(save.path,'postQC.rds'))
se.list <- lapply(se.list,function(se){
  se <- SCTransform(se,vst.flavor = "v2",assay = 'Spatial') 
  return(se)
})

# Batch correction(Harmony)---------------------------------------------------------------
var.features <- SelectIntegrationFeatures(object.list = se.list, nfeatures = 3000)
se <- merge(x = se.list[[1]], y = se.list[2:length(se.list)], merge.data = T)
VariableFeatures(se) <- var.features

# PCA
se <- RunPCA(se, assay = 'SCT', verbose = FALSE) 

# UMAP(no batch correction)
preBatch <- RunUMAP(se,reduction = "pca", dims = 1:30)

# Harmony batch correction
se <- RunHarmony(se,group.by.vars = "orig.ident", assay.use = "SCT", dims.use = 1:30, reduction = 'pca')
se <- RunUMAP(se,reduction = "harmony", dims = 1:30)

# save
saveRDS(preBatch, paste0(save.path,'preBatch.rds'))

# clustering---------------------------------------------------------------
se <- FindNeighbors(se, reduction = "harmony", dims = 1:30)
se <- FindClusters(se, resolution = seq(0.1,2,0.1)) 

# resolution = 0.8
Idents(se) <- se$seurat_clusters <- se$SCT_snn_res.0.8
se@misc$resolution <- 0.8

# save
saveRDS(se, paste0(save.path,'postClust.rds'))

# DE analysis---------------------------------------------------------------
se <- PrepSCTFindMarkers(se)
Idents(se) <- se$seurat_clusters
DEG <- FindAllMarkers(se,assay = 'SCT',logfc.threshold = 0.25, densify = T)