#5 Assortativity analysis


# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','SPATA2','semla'))

# Create STutility object---------------------------------------------------------------
se <- readRDS(paste0(save.path,'postAnnot.rds'))

ST_list <- lapply(c('Young1','Elderly1','Elderly2'),function(sample){
  sub_se <- subset(se,orig.ident==sample)
  infoTable <- data.frame(samples = paste0(sample,"/filtered_feature_bc_matrix.h5"),
                          spotfiles = paste0(sample,"/spatial/tissue_positions_list.csv"),
                          imgs = paste0(sample,"/spatial/spatial/tissue_hires_image.png"),
                          json = paste0(sample,"/spatial/spatial/scalefactors_json.json"))
  ST_obj <- ReadVisiumData(infoTable = infoTable)
  sub_se <- RenameCells(sub_se,new.names = gsub(colnames(sub_se), pattern = paste0('^',sample,'_'), replacement = ''))
  ST_obj <- SubsetSTData(ST_obj, spots = colnames(sub_se))
  all.equal(colnames(ST_obj), colnames(sub_se))
  ST_obj@meta.data <- sub_se@meta.data
  ST_obj <- LoadImages(ST_obj, time.resolve = F, verbose = F)
  MapLabels(ST_obj,column_name = "annotation")
  ST_obj <- SetIdent(ST_obj, value = "annotation")
  return(ST_obj)
})
names(ST_list) <- c('Young1','Elderly1','Elderly2')

# save
saveRDS(ST_list, paste0(save.path,'ST_list.rds'))


# Label Assortativity Analysis---------------------------------------------------------------

spata.list <- lapply(c('Young1','Elderly1','Elderly2'),function(sample){
  loadSpataObject(directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))
})
names(spata.list) <- c('Young1','Elderly1','Elderly2')
ROI_df <- lapply(c('Young1','Elderly1','Elderly2'),function(sample){
  lapply(c('ROI1','ROI2','ROI3','ROI4'),function(id){
    tmp <- getTrajectory(spata.list[[sample]],id = id)
    tmp@projection %>% 
      mutate(ROI = id) %>% 
      dplyr::select(barcodes, Sample = sample, ROI)
  }) %>% do.call(rbind,.)
}) %>% do.call(rbind,.) 

assort_res <- lapply(c('Young1','Elderly1','Elderly2'),function(sample){
  assort_res <- lapply(c('ROI1','ROI2','ROI3','ROI4'),function(id){
    ST_obj <- ST_list[[sample]]
    select_spot <- ROI_df %>% filter(Sample == sample & ROI == id) %>% pull(barcodes)
    res_lat <- ST_obj %>% 
      SubsetSTData(spots = select_spot) %>% 
      # SubsetSTData(idents = c('ZG','ZF','ZR')) %>% 
      RunLabelAssortativityTest(., column_name = "annotation")
    res_lat$ROI <- id
    return(res_lat)
  }) %>% do.call(rbind,.)
  assort_res$Sample <- sample
  return(assort_res)
}) %>% do.call(rbind,.)

## wilcox test
test_res <- lapply(c('ZG','ZF','ZR'), function(layer){
  young <- assort_res %>% filter(label == layer & Sample == 'Young1') %>% pull(avg_k_scaled)
  elderly <- assort_res %>% filter(label == layer & Sample != 'Young1') %>% pull(avg_k_scaled)
  res <- wilcox.test(young, elderly) 
  data.frame(Young_mean = mean(young, na.rm = T), Elderly_mean = mean(elderly,na.rm = T), P = res$p.value, Layer = layer)
}) %>% do.call(rbind,.)

# save
saveRDS(assort_res,paste0(save.path,'assort_res.rds'))