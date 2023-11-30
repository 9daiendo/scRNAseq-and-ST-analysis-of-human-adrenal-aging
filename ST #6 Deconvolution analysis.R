#6 Deconvolution analysis

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','spacexr','escape','msigdbr'))
Sys.setenv("OPENBLAS_NUM_THREADS"=8)

# Deconvolution analysis---------------------------------------------------------------
# prepare reference data (scRNAseq data)
SCdata <- readRDS('scRNAseq_postannotation.rds')
SCdata$cell_type <- SCdata$annotation

# prepare STdata
STdata <- readRDS(paste0(save.path,'postAnnot.rds'))

# Deconvolution
deconv_se <- list(0,0,0) ; names(deconv_se) <- unique(STdata$orig.ident)
for(sample in names(deconv_se)){
  # create spacexr object
  subSTdata <- subset(STdata, orig.ident == sample)
  
  subSCdata <- subset(SCdata, orig.ident == sample)
  SCreference <- Reference(counts = GetAssayData(subSCdata,assay = "RNA", slot = "counts"),
                           cell_types = subSCdata$cell_type)
  
  vis_coords <- GetTissueCoordinates(subSTdata,image = sample)
  VisiumData <- SpatialRNA(coords = vis_coords, counts = GetAssayData(subSTdata,slot = "counts", assay = "Spatial"))
  barcodes <- colnames(VisiumData@counts)

  # run RCTD
  myRCTD <- create.RCTD(VisiumData, SCreference, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  
  # Create variables from the myRCTD object to plot results
  barcodes <- colnames(myRCTD@spatialRNA@counts) # list of spatial barcodes
  weights <- myRCTD@results$weights # Weights for each cell type per barcode
  
  # Normalize per spot weights so cell type probabilities sum to 1 for each spot
  norm_weights <- normalize_weights(weights) 
  cell_type_names <- colnames(norm_weights) # List of cell types
  
  # add result to seurat object
  subSTdata@meta.data <- merge(subSTdata@meta.data,norm_weights,by=0,all=T) %>% column_to_rownames(var="Row.names")
  deconv_se[[sample]] <- subSTdata
}

# add result to seurat object
deconv_df <- lapply(names(deconv_se),function(x){
  tmp <- deconv_se[[x]]@meta.data[,cell_type_names]
  colnames(tmp) <- paste0('Prob_',colnames(tmp))
  return(tmp)
}) %>% do.call(rbind,.)
deconv_df <- deconv_df[Cells(STdata),]
df_celltype <- deconv_df %>% rownames_to_column(var = 'Cell') %>%
  pivot_longer(cols = -Cell, values_to = 'Prob', names_to = 'CellType') %>% 
  group_by(Cell) %>% arrange(-Prob) %>% dplyr::slice(1:1) %>% ungroup() %>%
  mutate(RCTD_CellType = gsub(CellType, pattern = 'Prob_', replacement = '')) %>%
  column_to_rownames(var = 'Cell')
deconv_df$RCTD_CellType <- df_celltype[rownames(deconv_df),'RCTD_CellType']
all.equal(rownames(deconv_df),Cells(STdata))
STdata@meta.data <- cbind(STdata@meta.data, deconv_df)

# save
saveRDS(STdata,paste0(save.path,"Deconv.rds"))


# ssGSEA---------------------------------------------------------------
# Gene Set
human.genes <- msigdbr(species = "Homo sapiens")
## HALLMARK geneset
hallmark <- human.genes %>% filter(gs_cat %in% c('H')) %>% split(x = .$gene_symbol, f = .$gs_name)
## Senescence geneset
senmayo <- read_tsv('SAUL_SEN_MAYO.v2023.1.Hs.tsv') %>%
  filter(STANDARD_NAME=='GENE_SYMBOLS') %>% pull(SAUL_SEN_MAYO) %>% strsplit(.,',')
names(senmayo) <- 'SenMayo'
## PKA
gobp <- human.genes %>% filter(gs_subcat %in% c('GO:BP')) %>% split(x = .$gene_symbol, f = .$gs_name)
pka <- gobp['GOBP_PROTEIN_KINASE_A_SIGNALING']

pathways.interest <- c(hallmark, senmayo, pka)

# run ssGSEA 
mat <- GetAssayData(se,slot = "counts", assay = "Spatial") # raw countを入力とする
ES <- enrichIt(obj = mat, gene.sets = pathways.interest, ssGSEA.norm = T, cores = 6,  min.size = 1)
se[['ssgsea']] <- CreateAssayObject(t(ES)[,Cells(se)])
DefaultAssay(se) <- 'ssgsea'

# save
saveRDS(se,paste0(save.path,"ssGSEA.rds"))