#4 select ROIs

# Setup---------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('Seurat','tidyverse','SPATA2'))

# create SPATA object--------------------------------------------------------
se <- readRDS(paste0(save.path,'postAnnot.rds'))

for(sample in c('Young1', 'Elderly1','Elderly2')){
  ## SPATA object
  spata_obj <- initiateSpataObject_10X(
    directory_10X = sample, # the directory from which to load the data
    sample_name = sample
  )
  ## metadata
  metadata <- se@meta.data %>% 
    filter(orig.ident == sample) %>% 
    rownames_to_column(var = 'barcodes') %>% 
    mutate(barcodes = gsub(barcodes, pattern = paste0("^",sample,'_'), replacement = ''))
  ## subset SPATA object
  spata_obj <- subsetByBarcodes(spata_obj, barcodes = metadata$barcodes)
  ## add metadata
  if(all.equal(spata_obj@fdata[[sample]]$barcodes, metadata$barcodes)){print("OK!!")}else{stop()}
  spata_obj <- addFeatures(object = spata_obj, feature_df = metadata, key_variable = 'barcodes', overwrite = T)
  ## denoising
  spata_obj <- runAutoencoderDenoising(spata_obj, activation = "selu", bottleneck = 56)
  # getActiveMatrixName(object = spata_obj)
  spata_obj <- setActiveMatrix(object = spata_obj, mtr_name = "denoised")
  # spata_obj <- setActiveMatrix(object = spata_obj, mtr_name = "scaled")
  # save spata object
  saveSpataObject(spata_obj, directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))
  }

# select ROIs--------------------------------------------------------
width <- '0.5mm'
## Young1------------
sample <- 'Young1'
spata_obj <- loadSpataObject(directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))
### ROI1
id <- 'ROI1'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "1mm", y = "2.8mm"),
                                  end = c(x = "2.1mm", y = "2.2mm"))
### ROI2
id <- 'ROI2'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "2.4mm", y = "4.6mm"),
                                  end = c(x = "3.3mm", y = "3.3mm"))
### ROI3
id <- 'ROI3'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "3.8mm", y = "6mm"),
                                  end = c(x = "4.6mm", y = "5mm"))
### ROI4
id <- 'ROI4'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "6.2mm", y = "7.1mm"),
                                  end = c(x = "6.7mm", y = "5.6mm"))

### save
saveSpataObject(spata_obj, directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))


## Elderly1------------
sample <- 'Elderly1'
spata_obj <- loadSpataObject(directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))
### ROI1
id <- 'ROI1'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "1.8mm", y = "3.3mm"),
                                  end = c(x = "2.8mm", y = "3.7mm"))
### ROI2
id <- 'ROI2'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "2.5mm", y = "2.2mm"),
                                  end = c(x = "3.3mm", y = "3mm"))
### ROI3
id <- 'ROI3'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "3.7mm", y = "1.2mm"),
                                  end = c(x = "4.3mm", y = "2.1mm"))
### ROI4
id <- 'ROI4'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "4.9mm", y = "0.8mm"),
                                  end = c(x = "5.2mm", y = "1.7mm"))

### save
saveSpataObject(spata_obj, directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))


## Elderly2------------
sample <- 'Elderly2'
spata_obj <- loadSpataObject(directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))
### ROI1
id <- 'ROI1'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "1.8mm", y = "2.8mm"),
                                  end = c(x = "2.6mm", y = "2mm"))
### ROI2
id <- 'ROI2'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "2.9mm", y = "4.1mm"),
                                  end = c(x = "3.5mm", y = "3.7mm"))
### ROI3
id <- 'ROI3'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "3.8mm", y = "5.4mm"),
                                  end = c(x = "4.5mm", y = "4.8mm"))
### ROI4
id <- 'ROI4'
spata_obj <- addSpatialTrajectory(spata_obj, id = id, overwrite = T,
                                  width = width,
                                  start = c(x = "5.1mm", y = "6.7mm"),
                                  end = c(x = "5.8mm", y = "6.1mm"))

### save
saveSpataObject(spata_obj, directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))

