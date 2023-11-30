#6 Trajectory analysis

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','SeuratWrappers','SeuratDisk',
                          'tidyverse','destiny','SingleCellExperiment','slingshot','org.Hs.eg.db','escape','msigdbr'))


# subset Young adrenocortical cells--------------------------------------------------
se <- readRDS(paste0(save.path,'postAnnot.rds'))
se <- subset(se, Age == 'Young' & annotation %in% c('ZG','ZF','ZR'))
Idents(se) <- se$annotation

# Diffusion map--------------------------------------------------
DefaultAssay(se) <- 'SCT'
se <- se %>% RunPCA(assay = 'SCT')
dm <- DiffusionMap(Embeddings(se, reduction = "pca")[,1:30], n_eigs = 2) 
se[["dc"]] <- CreateDimReducObject(embeddings = apply(as.matrix(as.data.frame(dm)[,1:2]),MARGIN = 2, rescale),
                                   key = "DC_", assay = 'RNA')

# Slingshot--------------------------------------------------
sce <- SingleCellExperiment(list(counts = se@assays$RNA@counts,
                                 normcounts = se@assays$RNA@data))
reducedDims(sce) <- list(pca = Embeddings(se, reduction = "pca"), 
                         umap = Embeddings(se, reduction = "umap"),
                         dc = Embeddings(se, reduction = "dc")
)
colData(sce)$annotation <- se$annotation

# Run slingshot
pto <- slingshot(sce, clusterLabels = sce$annotation, reducedDim = 'dc', start.clus = 'ZG')
sds <- as.SlingshotDataSet(pto)
se$Pseudotime <- colData(pto)$slingPseudotime_1

# save
saveRDS(se,paste0(save.path,'TrajectoryYoungCortex.rds'))
saveRDS(pto,paste0(save.path,'TrajectoryObj.rds'))

# ssGSEA---------------------------------------
# Gene Set
human.genes <- msigdbr(species = "Homo sapiens")
## HALLMARK geneset
hallmark <- human.genes %>% filter(gs_cat %in% c('H')) %>% split(x = .$gene_symbol, f = .$gs_name)
## Senescence geneset
senmayo <- read_tsv('Rawdata/SAUL_SEN_MAYO.v2023.1.Hs.tsv') %>%
  filter(STANDARD_NAME=='GENE_SYMBOLS') %>% pull(SAUL_SEN_MAYO) %>% strsplit(.,',')
names(senmayo) <- 'SenMayo'
other_sen <- human.genes %>% 
  filter(gs_cat %in% c('C2')) %>% 
  filter(gs_name %in% c(
    'FRIDMAN_SENESCENCE_UP',
    'REACTOME_CELLULAR_SENESCENCE',
    'REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP'
  )) %>% 
  split(x = .$gene_symbol, f = .$gs_name)
## PKA
gobp <- human.genes %>% filter(gs_subcat %in% c('GO:BP')) %>% split(x = .$gene_symbol, f = .$gs_name)
pka <- gobp['GOBP_PROTEIN_KINASE_A_SIGNALING']

pathways.interest <- c(hallmark, senmayo, pka)

# run ssGSEA 
mat <- GetAssayData(se,slot = "counts", assay = "RNA") # input raw count
ES <- enrichIt(obj = mat, gene.sets = pathways.interest, ssGSEA.norm = T, cores = 6,  min.size = 1)
se[['ssgsea']] <- CreateAssayObject(t(ES)[,Cells(se)])
DefaultAssay(se) <- 'ssgsea'

# save
saveRDS(se,paste0(save.path,'ssGSEA_Young.rds'))


# convert h5ad---------------------------------------
DefaultAssay(se) <- 'RNA'
Idents(se) <- se$annotation

# organize metadata
se@meta.data <- se@meta.data %>% 
  dplyr::select(orig.ident, nCount_RNA, nFeature_RNA,percent.mt,S.Score,G2M.Score,Phase,Cycle.Score,Age,annotation,sub_annotation,Pseudotime) %>% 
  mutate(Age = as.character(Age),annotation = as.character(annotation),sub_annotation = as.character(sub_annotation))
se$barcodes <- Cells(se)

# organize dim reduction
se[['pca']] <- CreateDimReducObject(embeddings = se@reductions$pca@cell.embeddings, assay = 'RNA')
se[['umap']] <- CreateDimReducObject(embeddings = se@reductions$umap@cell.embeddings, assay = 'RNA')
se@reductions

# narrow to top3000 var genes
all.data <- readRDS(paste0(save.path,'postAnnot.rds'))
se <- se[all.data@assays$SCT@var.features,]
dim(se)

# delete SCT assay
se[['SCT']] <- NULL

# convert to h5ad
SaveH5Seurat(se, filename = paste0(save.path,'Young_cortex_3000genes.h5Seurat'), overwrite = T)
Convert(paste0(save.path,'Young_cortex_3000genes.h5Seurat'), dest = "h5ad", overwrite = T, assay = 'RNA')



