#2 Merge data and clustering

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','harmony'))

# merge objects and rename cells---------------------------------------------------------------
se.list <- readRDS(paste0(save.path,'postQC.rds'))
se.list <- lapply(se.list,function(se){
  se <- RenameCells(object = se, add.cell.id = unique(se$orig.ident))
})

# Batch correction(Harmony)---------------------------------------------------------------
# integration
se.list <- lapply(se.list,function(se){
  # Log Normalize(need before cellcyclescoring)
  se <- NormalizeData(se)
  # Cel cycle score
  se <- CellCycleScoring(se, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)
  se$Cycle.Score <- se$S.Score - se$G2M.Score
  # SCTransform
  se <- SCTransform(se,vst.flavor = "v2",vars.to.regress = c("Cycle.Score", "percent.mt")) 
  return(se)
})
var.features <- SelectIntegrationFeatures(object.list = se.list, nfeatures = 3000)
se <- merge(x = se.list[[1]], y = se.list[2:length(se.list)], merge.data = T)
VariableFeatures(se) <- var.features
se$Age <- factor(ifelse(se$orig.ident == 'Young1', 'Young', 'Elderly'),levels = c('Elderly','Young'))

# PCA
se <- RunPCA(se, assay = 'SCT', verbose = FALSE) 
ElbowPlot(se,ndims = 50)

# UMAP(no batch correction)
preBatch <- RunUMAP(se,reduction = "pca", dims = 1:30)

# Harmony batch correction
se <- RunHarmony(se,group.by.vars = "orig.ident", assay.use = "SCT", dims.use = 1:30, reduction = 'pca')
se <- RunUMAP(se,reduction = "harmony", dims = 1:30)

# save
saveRDS(preBatch, paste0(save.path,'preBatch.rds'))

# Clustering---------------------------------------------------------------
se <- FindNeighbors(se, reduction = "harmony", dims = 1:30)
se <- FindClusters(se, resolution = seq(0.1,1.2,0.1)) 

# resolution = 0.4
Idents(se) <- se$seurat_clusters <- se$SCT_snn_res.0.4

# save
saveRDS(se, paste0(save.path,'postClust.rds'))

# DE analysis---------------------------------------------------------------
DefaultAssay(se) <- 'SCT'
se <- PrepSCTFindMarkers(se)

# DE analysis
Idents(se) <- se$seurat_clusters
DEG <- FindAllMarkers(se,assay = 'SCT',logfc.threshold = 0.25, densify = T)
