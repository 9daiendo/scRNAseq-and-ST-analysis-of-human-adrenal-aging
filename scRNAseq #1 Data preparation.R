#1 Data preparation

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','patchwork','celda','singleCellTK','scDblFinder','SeuratWrappers','reticulate'))

# remove ambient RNA(decontX)---------------------------------------------------------------
se.list <- lapply(samples,function(sample){
  # Create a SingleCellExperiment object and run decontX
  counts <- Read10X_h5(filename = paste0(sample,'filtered_feature_bc_matrix.h5'))
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- celda::decontX(sce)
  
  # Create a Seurat object from a SCE with decontX results
  se <- CreateSeuratObject(round(celda::decontXcounts(sce)), project = sample)
  return(se)
})

# Doublet detection(scDblfinder)---------------------------------------------------------------
se.list <- lapply(se.list,function(se){
  counts <- GetAssayData(se,assay = 'RNA', slot = 'count')
  sce <- SingleCellExperiment(list(counts = counts))
  ## scDblfinder
  sce <- scDblFinder(sce, BPPARAM = BiocParallel::MulticoreParam(6,RNGseed=1234))
  se$doublet_class <- sce$scDblFinder.class
  se$doublet_score <- sce$scDblFinder.score
  se$percent.mt <- PercentageFeatureSet(se,pattern = "^MT-") # mitochondrial gene fraction
  return(se)
})

# save
saveRDS(se.list, paste0(save.path,'preQC.rds'))

# Filtering---------------------------------------------------------------
feature.max <- 8000 ; feature.min <- 200 ; mt.cut <- 20

QCtable <- lapply(unique(plot.df$orig.ident),function(x){
  tmp <- plot.df %>% filter(orig.ident==x)
  data.frame(
    Sample = x,
    Total = nrow(tmp),
    Doublet = sum(tmp$doublet_class=='doublet'),
    FeatureHigh = sum(tmp$nFeature_RNA >= feature.max),
    FeatureLow = sum(tmp$nFeature_RNA <= feature.min),
    MitoHigh = sum(tmp$percent.mt >= mt.cut),
    AllOK = sum(tmp$doublet_class=='singlet' & 
                  tmp$nFeature_RNA < feature.max &
                  tmp$nFeature_RNA > feature.min &
                  tmp$percent.mt < mt.cut))
}) %>% bind_rows()

se.list <- lapply(se.list,function(se){
  se <- subset(se,
               nFeature_RNA > feature.min & 
               nFeature_RNA < feature.max & 
               percent.mt < mt.cut &
               doublet_class == "singlet")
  return(se)
})

# save
saveRDS(se.list, paste0(save.path,'postQC.rds'))
