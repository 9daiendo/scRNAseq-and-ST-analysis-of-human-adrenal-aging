# GTEx analysis

# setup--------------------------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
set.seed(1234)

wdir <- 'set your working directory'
setwd(wdir)
save.path <- 'object folder path'

# load library
easypackages::libraries(c('recount3','edgeR','tidyverse','org.Hs.eg.db','patchwork','fgsea','msigdbr'))

# Data preparation--------------------------------------------------------------------------------
# Get GTEx data 
df <- recount3::create_rse_manual(
  project = "ADRENAL_GLAND",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene"
)

# Sample QC
## age group
df@colData$age <- factor(ifelse(df@colData$gtex.age%in%c("60-69","70-79"),"Elderly","Young"),levels = c("Elderly","Young"))
## sex group
df@colData$sex <- factor(ifelse(df@colData$gtex.sex==1,"Male","Female"), levels = c("Male","Female"))
## gtex metadata
df@colData[,grep('gtex.',colnames(df@colData),value = T)]
## Exclude low RNA quality sample
df <- df[,df$gtex.smafrze=='RNASEQ']
df <- df[,df$gtex.smrin>=6]

# Modify gene expression matrix
## Omit duplicate genes
gene_ord <- rowSums(assay(df,'raw_counts')) %>% sort(decreasing = T) %>% names()
df <- df[gene_ord,]
gene_df <- rowRanges(df) %>% as.data.frame() %>% 
  distinct(gene_name, .keep_all = T)
df <- df[gene_df$gene_id,]
rownames(df) <- gene_df$gene_name
## Scaling (According to the manual of recount3)
assay(df, "counts") <- transform_counts(df)

# Sample labeling
df@colData$Label[df@colData$age == 'Elderly'] <- paste0('Elderly',base::seq(1,sum(df@colData$age == 'Elderly'),1))
df@colData$Label[df@colData$age == 'Young'] <- paste0('Young',seq(1,sum(df@colData$age == 'Young'),1))
colnames(df) <- df@colData$Label

# Create edgeR object
count_df <- assay(df,'counts')
gene_df <- rowRanges(df) %>% as.data.frame()
sample_df <- as.data.frame(df@colData)
## Check the sample order
all.equal(rownames(count_df),rownames(gene_df))
all.equal(colnames(count_df),rownames(sample_df))

d <- DGEList(counts = count_df, genes = gene_df, samples = sample_df)

# Sample info
d$samples %>% 
  dplyr::select(age,sex,gtex.age,gtex.smrin) %>% 
  tableone::CreateTableOne(data = ., strata = "age", argsNonNormal = 'gtex.smrin')


# QC, Normalization--------------------------------------------------------------------------------
# Exclude low expression genes
keep.exprs <- filterByExpr(d)
table(keep.exprs)
## Check steroidogenic genes
steroidogenic <- c("CYP11A1","CYP11B1","CYP11B2","CYP17A1","CYP21A2","HSD3B2","CYB5A","SULT2A1")
steroidogenic[!steroidogenic%in%names(keep.exprs)[keep.exprs]]

d <- d[keep.exprs,,keep.lib.sizes = F]
dim(d) 

# TMM Normalising gene expression distributions
d <- calcNormFactors(d, method = "TMM")
d$samples$norm.factors 
cpm <- edgeR::cpm(d, log=F) # CPM
lcpm <- edgeR::cpm(d, log=T) # log2CPM

# save objects
saveRDS(d, paste0(save.path,'edgeRobj.rds'))

# DE analysis(EdgeR)--------------------------------------------------------------------------------
# setting comparison
## elderly vs young
group1 <- 'Elderly' ; group2 <- 'Young'
contrast <- paste0(group1,'-',group2)
select <- d$samples %>% 
  filter(age %in% c(group1,group2)) %>% rownames()
group <- factor(d$samples[select,'age'],levels = c(group2,group1))

## elderly male vs young female
group1 <- 'ElderlyMale' ; group2 <- 'YoungFemale'
contrast <- paste0(group1,'-',group2)
group1_select <- d$samples %>% 
  filter(age == 'Elderly' & sex == 'Male') %>% rownames()
names(group1_select) <- rep(group1,length(group1_select))
group2_select <- d$samples %>% 
  filter(age == 'Young' & sex == 'Female') %>% rownames()
names(group2_select) <- rep(group2,length(group2_select))

# DE analysis
select <- c(group1_select,group2_select)
group <- factor(names(select),levels = c(group2,group1))
design <- model.matrix(~group)
select_d <- d[,select] ; select_d <- estimateDisp(select_d,design)
select_lcpm <- lcpm[,select]

fit <- glmFit(select_d, design)
lrt <- glmLRT(fit,coef = 2) 
result <- topTags(lrt, n = Inf, sort.by="logFC", adjust.method="BH")
result <- result %>% as.data.frame() %>% 
  dplyr::select(gene_name, log2FC = logFC, p_val_adj = FDR, everything())

# save object
saveRDS(result, paste0(save.path,'DEA_res.rds'))

# view result
result %>% 
  filter(abs(log2FC) > 0.25 & p_val_adj < 0.05) %>%
  arrange(-log2FC) %>% 
  view()

# GSEA---------------------------------------------------------------
# Gene Set
human.genes <- msigdbr(species = "Homo sapiens")
## HALLMARK geneset
hallmark <- human.genes %>% 
  filter(gs_cat %in% c('H')) %>% 
  split(x = .$gene_symbol, f = .$gs_name)
## Senescence geneset
senmayo <- read_tsv('SAUL_SEN_MAYO.v2023.1.Hs.tsv') %>% # download from https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SAUL_SEN_MAYO
  filter(STANDARD_NAME=='GENE_SYMBOLS') %>% 
  pull(SAUL_SEN_MAYO) %>% strsplit(.,',')
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
gobp <- human.genes %>% 
  filter(gs_subcat %in% c('GO:BP')) %>% 
  split(x = .$gene_symbol, f = .$gs_name)
pka <- gobp['GOBP_PROTEIN_KINASE_A_SIGNALING']

pathways.interest <- c(hallmark,pka)
pathways.interest <- c(senmayo,other_sen)

# GSEA
stats <- result %>% 
  arrange(desc(log2FC)) %>% 
  dplyr::select(gene_name, log2FC) %>% 
  deframe()
gsea_res <- fgsea(pathways = pathways.interest, stats = stats, nperm = 10000) 

# save
saveRDS(gsea_res, paste0(save.path,'GSEA_res.rds'))

