#8 Ligand_receptor analysis


# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c("tidyverse","Seurat","liana","CrossTalkeR"))

# run liana------------------------------------------------------------------------------------
se <- readRDS(paste0(save.path,'postAnnot.rds'))

for(i in c('Young','Elderly')){
  sub_se <- subset(se,Age == i)
  Idents(sub_se) <- sub_se$annotation
  DefaultAssay(sub_se) <- 'RNA'
  liana_results <- sub_se %>% liana_wrap() %>% liana_aggregate()
  write.csv(liana_results, file = paste0(save.path,'LigandReceptor_',i,'.csv'))
}

# CrossTalkeR------------------------------------------------------------------------------------
dir.create(paste0(save.path,"CrossTalkeR")) # directory for crosstalker

for(i in c('Young','Elderly')){
  liana_results <- read.csv(paste0(save.path,'LigandReceptor_',i,'.csv'))

  liana_results <- liana_results %>% 
    filter(source %in% c('ZG','ZF','ZR','Macrophage')) %>% 
    filter(target %in% c('ZG','ZF','ZR','Macrophage')) 
    
  df <- liana_results %>% 
    dplyr::rename(ligand=ligand.complex, receptor=receptor.complex) %>% 
    mutate(id = fct_inorder(paste0(ligand, " -> ", receptor))) %>% 
    filter(cellphonedb.pvalue < 0.05) %>% # cellphonedb pvalue
    mutate(type_gene_A = "Ligand", type_gene_B = "Receptor") %>% 
    rename(gene_A = ligand, gene_B = receptor, MeanLR = sca.LRscore) %>% 
    dplyr::arrange(source, target, aggregate_rank) %>% 
    dplyr::select(source, target, gene_A, gene_B, type_gene_A, type_gene_B, MeanLR)
  write.csv(df, file = paste0(save.path,'CrossTalkeR/',i,'_df.csv'))
}

paths <- c('CTR' = paste0(save.path,'CrossTalkeR/Young_df.csv'), 
           'EXP' = paste0(save.path,'CrossTalkeR/Elderly_df.csv'))
data <- generate_report(paths,
                        out_path = paste0(save.path,"CrossTalkeR/"), # 絶対パス
                        threshold = 0,
                        out_file = 'Result.html',
                        output_fmt = "html_document",
                        report = F)
