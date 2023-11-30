#5 SCENIC

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','SCENIC','AUCell'))

# load pySCENIC data---------------------------------------
load(file = paste0(save.path,"SCENIC.rda"))
se <- readRDS(paste0(save.path,'postAnnot.rds'))
se <- subset(se,annotation%in%c('ZG','ZF','ZR'))

# add AUC assay
AUCmat <- AUCell::getAUC(regulonAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))
names(regulons) <- gsub("[(+)]", "", names(regulons))
se[['AUC']] <- CreateAssayObject(data = AUCmat)

# regulon data
regulon_df <- lapply(1:length(regulons),function(i){
  tibble(TF = names(regulons[i]), regulon = regulons[[i]])
}) %>% bind_rows()


# Cell-type specific regulators (RSS)---------------------------------------
se$Age <- ifelse(se$orig.ident == "Young1", "Young", "Elderly")
Idents(se) <- paste0(se$Age,'_',se$annotation) %>% 
  factor(.,levels = c(paste0('Young_',c('ZG','ZF','ZR')),paste0('Elderly_',c('ZG','ZF','ZR'))))

rss <- calcRSS(AUC = se@assays$AUC@data, 
               cellAnnotation = Idents(se)
               ) %>% 
  as_tibble(rownames = "TF") %>% 
  pivot_longer(cols = -TF, names_to = "annotation", values_to = "RSS") %>% 
  drop_na(RSS) %>% group_by(annotation) %>% mutate(rank = rank(-RSS)) %>% ungroup()

# save
saveRDS(rss,paste0(save.path,'RSS_df.rds'))