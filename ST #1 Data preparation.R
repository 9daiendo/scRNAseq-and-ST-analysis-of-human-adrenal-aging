#1 Data preparation, Quality control

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse'))

# make seurat object------------------------------------------------------------
se.list <- lapply(c('Young1', 'Elderly1','Elderly2'),function(sample){
  se <- Load10X_Spatial(data.dir = sample,
                        filename = "filtered_feature_bc_matrix.h5",
                        slice = sample,
                        image = Read10X_Image(image.dir = paste0(sample,"/spatial")))
  annot_df <- read.csv(paste0(sample,"/Histological_annotation.csv"), row.names = 1)
  se$Histological_annotation <- annot_df[colnames(se),]
  se$orig.ident <- sample
  se$percent.mt <- PercentageFeatureSet(se, pattern = "^MT-")
  se <- CellCycleScoring(object = se, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)
  se <- RenameCells(se,add.cell.id = sample)
  return(se)
})

# save
saveRDS(se.list, paste0(save.path,'preQC.rds'))

# Filtering------------------------------------------------------------
feature.min <- 500

QCtable <- lapply(unique(plot.df$orig.ident),function(x){
  tmp <- plot.df %>% filter(orig.ident==x)
  data.frame(
    Sample = x,
    Total = nrow(tmp),
    FeatureLow = sum(tmp$nFeature_Spatial <= feature.min),
    HistExclude = sum(tmp$Histological_annotation == 'Exclude'),
    AllOK = sum(tmp$nFeature_Spatial > feature.min & 
                  tmp$Histological_annotation != 'Exclude'))
}) %>% bind_rows()

# Filtering
se.list <- lapply(se.list,function(x){subset(x,nFeature_Spatial>feature.min & Histological_annotation != "Exclude")})

# save
saveRDS(se.list, paste0(save.path,'postQC.rds'))
