#4 Compare Young vs Elderly

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','fgsea','msigdbr','escape','org.Hs.eg.db'))

# compare cell composition---------------------------------------------------------------
se <- readRDS(paste0(save.path,'postAnnot.rds'))

table(se@meta.data$Age)
se@meta.data %>% as.data.frame() %>% 
  xtabs(~Age + annotation, data = .) 

tab <- se@meta.data %>% as.data.frame() %>% 
  filter(annotation %in% c('ZG','ZF','ZR')) %>% 
  mutate(annotation = factor(annotation,levels = c('ZG','ZF','ZR'))) %>% 
  mutate(Age = factor(Age, levels = c('Young', 'Elderly'))) %>% 
  xtabs(~Age + annotation, data = .) 
fisher.test(tab) 
RVAideMemoire::fisher.multcomp(tab, p.method = "BH") 

# DE analysis Young vs Elderly---------------------------------------------------------------
Idents(se) <- se$Age
# all
DEG <- FindMarkers(se,assay = 'SCT',
                   # test.use = 'MAST',latent.vars = 'orig.ident', # batch correct
                   ident.1 = 'Elderly',ident.2 = 'Young',
                   logfc.threshold = 0.25,densify = T)
DEG %>% 
  filter(p_val_adj < 0.05) %>% 
  arrange(-avg_log2FC) %>% 
  rownames_to_column(var = 'gene') %>% 
  xlsx::write.xlsx(file = paste0(save.path,'DEG_ElderlyvsYoung.xlsx'),row.names = F)

# each annotation cluster
for(i in levels(se$annotation)){
  sub_se <- subset(se,annotation==i)
  sub_se <- PrepSCTFindMarkers(sub_se)
  DEG <- FindMarkers(sub_se, assay = 'SCT',
                     # test.use = 'MAST',latent.vars = 'orig.ident', recorrect_umi = F, # batch correct
                     ident.1 = 'Elderly',ident.2 = 'Young',
                     logfc.threshold = 0.25, densify = T)
  DEG %>% 
    filter(p_val_adj < 0.05) %>% arrange(-avg_log2FC) %>% rownames_to_column(var = 'gene') %>% 
    xlsx::write.xlsx(file = paste0(save.path,'DEG_ElderlyvsYoung_',i,'.xlsx'),row.names = F)
}


# GSEA(Hallmark)---------------------------------------------------------------
# Gene Set
human.genes <- msigdbr(species = "Homo sapiens")
## HALLMARK geneset
hallmark <- human.genes %>% filter(gs_cat %in% c('H')) %>% split(x = .$gene_symbol, f = .$gs_name)
## Senescence geneset
senmayo <- read_tsv('SAUL_SEN_MAYO.v2023.1.Hs.tsv') %>%
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


# GSEA
df <- lapply(c('ZG','ZF','ZR'),function(cluster){
  sub_se <- subset(se, annotation == cluster)
  # fold change
  fc_res <- AverageExpression(sub_se,assays = 'SCT', features = rownames(sub_se), 
                              return.seurat = F, group.by = 'Age', slot = 'data')[[1]] %>% 
    as.data.frame() %>% mutate(avg_log2FC = log2(Elderly+1)-log2(Young+1)) %>% 
    rownames_to_column(var = 'gene')
  # GSEA
  stats <- fc_res %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene, avg_log2FC) %>% deframe()
  gsea_res <- fgsea(pathways = pathways.interest, stats = stats, nperm = 10000) 
  gsea_res$cluster <- cluster
  return(gsea_res)
  }) %>% do.call(rbind,.)

# save
saveRDS(df, paste0(save.path,'GSEA_ElderlyvsYoung.rds'))