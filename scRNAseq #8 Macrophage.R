#9 Macrophage

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','fgsea','msigdbr','escape','org.Hs.eg.db'))

# GSEA---------------------------------------------------------------
se <- readRDS(paste0(save.path,'postAnnot.rds'))
se <- subset(se, annotation == 'Macrophage')

# Gene Set
human.genes <- msigdbr(species = "Homo sapiens")
## HALLMARK geneset
hallmark <- human.genes %>% filter(gs_cat %in% c('H')) %>% split(x = .$gene_symbol, f = .$gs_name)
## M1, M2 geneset
m1 <- cmapR::parse_gmt('Rawdata/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.v2023.1.Hs.gmt')
m1 <- m1[[1]]$entry[m1[[1]]$entry%in%rownames(se)]
m2 <- cmapR::parse_gmt('Rawdata/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.v2023.1.Hs.gmt')
m2 <- m2[[1]]$entry[m2[[1]]$entry%in%rownames(se)]
M1M2 <- list(M1_macrophage = m1, M2_macrophage = m2)

pathways.interest <- c(hallmark,M1M2)

# calculate fold change
DefaultAssay(se) <- 'SCT'
Idents(se) <- se$Age
fc_res <- AverageExpression(se,assays = 'SCT', features = rownames(se), group.by = 'Age',
                            return.seurat = F, slot = 'data')[[1]] %>% 
  as.data.frame() %>% 
  mutate(avg_log2FC = log2(Elderly+1)-log2(Young+1)) %>% 
  rownames_to_column(var = 'gene')

# GSEA
stats <- fc_res %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC) %>% deframe()
gsea_res <- fgsea(pathways = pathways.interest, stats = stats, nperm = 10000) 

# save
saveRDS(gsea_res, paste0(save.path,'GSEA_Mac.rds'))