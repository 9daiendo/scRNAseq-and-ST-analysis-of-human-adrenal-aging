# bulkRNAseq------------------------------------------------------------
## setup----
# load library
easypackages::libraries(c('recount3','edgeR','tidyverse','org.Hs.eg.db','patchwork','fgsea','msigdbr'))
# load objects
save.path <- 'object folder path'
d <- readRDS(paste0(save.path,'edgeRobj.rds'))
lcpm <- edgeR::cpm(d, log=T) 
result <- readRDS(paste0(save.path,'DEA_res.rds'))
gsea_res <- readRDS(paste0(save.path,'GSEA_res.rds'))
# color palette
age_pal <- c('#C77CFF','#7CAE00') ; names(age_pal) <- c('Elderly','Young')
sex_pal <- c('#00BFC4','#F8766D') ; names(sex_pal) <- c('Male','Female')
# figure save folder
fig.path <- 'figure folder path'

## Violin plots of steroidogenic genes----
exp_df <- lcpm %>% t() %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample") %>%
  mutate(group = factor(d$samples$age))

steroidogenic <- c("CYP11A1","CYP11B1","CYP11B2","CYP17A1","CYP21A2","HSD3B2","CYB5A","SULT2A1")
gene_list <- steroidogenic

for(gene in gene_list){
  log2FC <- result[gene,"log2FC"] ; p_val_adj <- result[gene,"p_val_adj"]
  lab <- ifelse(abs(log2FC) > 0.25 & p_val_adj < 0.05, "*","")
  plot_df <- exp_df %>% dplyr::select(gene = all_of(gene), group, sample)
  plot_df %>%  ggplot(aes(x = group, y = gene)) + 
    geom_violin(aes(fill = group), color = NA, alpha = 0.6) + 
    scale_fill_manual(values = age_pal, name = "age") +
    geom_boxplot(width = 0.07, outlier.colour = NA, position = position_dodge(width = 0.9)) +
    scale_x_discrete(limits = c('Young','Elderly')) +
    ggsignif::geom_signif(comparisons = list(levels(plot_df$group)), 
                          map_signif_level = F, 
                          annotations = lab, textsize = 8, fontface = "bold", vjust = 0.5) +
    ylim(c(NA,max(plot_df$gene)*1.1)) +
    theme_classic() + labs(title = gene, x = "", y = "Expression Level") +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          legend.position = "none",
          aspect.ratio = 1)
  ggsave(paste0(fig.path,"ViolinPlot_",gene,".pdf"),height = 4, width = 4, dpi = 300)
}

## GSEA bar plot----
plot.df <- gsea_res %>% 
  group_by(NES > 0) %>% 
  arrange(-abs(NES)) %>% 
  mutate(order = row_number()) %>% 
  ungroup() %>% 
  mutate(pathway = gsub(pathway,pattern = 'HALLMARK_',replacement = '')) %>% 
  mutate(pathway = gsub(pathway,pattern = 'GOBP_',replacement = '')) %>% 
  mutate(pathway = str_to_title(pathway)) %>% 
  mutate(pathway = gsub(pathway,pattern = '_',replacement = ' ')) %>% 
  filter(padj < 0.05) %>% group_by(NES>0) %>% 
  top_n(n = 15, wt = abs(NES)) 
plot.df %>% 
  ggplot(aes(x = NES, y = fct_reorder(pathway, NES))) + 
  geom_col(aes(fill = -log10(padj))) +
  scale_fill_viridis_c(option = 'viridis') +
  # xlim(c(-max(abs(plot.df$NES)),max(abs(plot.df$NES)))) +
  theme_bw() +
  scale_y_discrete(position = "right") +
  theme(aspect.ratio = 2, axis.text = element_text(color = 'black'),
        panel.grid =  element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())

ggsave(paste0(fig.path,'GSEA_barplot.pdf'), height = 5, width = 6, dpi = 300)

## GSEA plot----
select_pathway <- 'HALLMARK_WNT_BETA_CATENIN_SIGNALING'
plotEnrichment(pathways.interest[[select_pathway]], stats, gseaParam = 1, ticksSize = 0.3) +
  xlab('Rank') + ylab('Enrichment Score') + 
  theme(aspect.ratio = 0.8, 
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) 

ggsave(paste0(fig.path,'GSEA_',select_pathway,'.pdf'),height = 4,width = 5,dpi = 300)

# scRNAseq------------------------------------------------------------
## setup----
# load library
easypackages::libraries(c('Seurat','tidyverse','patchwork','harmony',
                          'ComplexHeatmap','circlize','RColorBrewer', 'cowplot'))

# load objects
save.path <- 'object folder path'
preQC <- readRDS(paste0(save.path,'preQC.rds'))
postQC <- readRDS(paste0(save.path,'postQC.rds'))
preBatch<- readRDS(paste0(save.path,'preBatch.rds'))
postClust <- readRDS(paste0(save.path,'postClust.rds'))
postAnnot <- readRDS(paste0(save.path,'postAnnot.rds'))
gsea_res <- readRDS(paste0(save.path,'GSEA_ElderlyvsYoung.rds'))
rss <- readRDS(paste0(save.path,'RSS_df.rds'))
Traj <- readRDS(paste0(save.path,'TrajectoryYoungCortex.rds'))
pto <- readRDS(paste0(save.path,'TrajectoryObj.rds'))
ssGSEA_Young <- readRDS(paste0(save.path,'ssGSEA_Young.rds'))
LRobj <- readRDS(paste0(save.path,'CrossTalkeR/LR_data_final.Rds'))
GSEA_Mac <- readRDS(paste0(save.path,'GSEA_Mac.rds'))

# color palette
age_pal <- c('#C77CFF','#7CAE00') ; names(age_pal) <- c('Elderly','Young')
annot_pal <- c('Capsule' = '#66C2A5', 'ZG' = '#FC8D62', 'ZF' = '#8DA0CB', 'ZR' = '#E78AC3',
               'Medulla' = '#A6D854', 'Endothelial' = '#FFD92F',
               'Fibroblast' = '#E5C494',"VSM" = '#B3B3B3',
               'Adipocyte' = '#B15928', "Macrophage" = '#FFFF99', "Tcell" = '#6A3D9A',
               "Bcell" = '#CAB2D6', 'Neutrophil' = '#FDBF6F'
)


# figure save folder
fig.path <- 'figure folder path'

## QC plot----
feature.max <- 8000 ; feature.min <- 200 ; mt.cut <- 20
se.list <- preQC
plot.df <- lapply(se.list,function(se){return(se@meta.data)}) %>% bind_rows()
plot.df$orig.ident <- factor(plot.df$orig.ident, levels = c('Young1','Elderly1','Elderly2'))
plot.set <- theme_classic() + 
  theme(aspect.ratio = 1, axis.title = element_blank(), 
        axis.text.x = element_text(size = 10), axis.ticks.x = element_blank(), 
        title = element_text(face = "bold"), legend.position = "right")
p1 <- plot.df %>% ggplot(aes(x = orig.ident, y = nCount_RNA)) +
  labs(title = "#UMIs") + geom_violin(color = "black", fill = "lightgrey") +
  plot.set
p2 <- plot.df %>% ggplot(aes(x = orig.ident, y = nFeature_RNA)) + 
  labs(title = "#Genes") + 
  geom_violin(color = "black", fill = "lightgrey") + 
  geom_hline(yintercept = feature.max, linetype = "dashed", color = "red") +
  geom_hline(yintercept = feature.min, linetype = "dashed", color = "red") +
  plot.set
p3 <- plot.df %>% ggplot(aes(x = orig.ident, y = percent.mt)) + 
  labs(title = "Percent.mt") + 
  geom_violin(color = "black", fill = "lightgrey") + 
  geom_hline(yintercept = mt.cut, linetype = "dashed", color = "red") +
  plot.set
wrap_plots(p1,p2,p3) + 
  patchwork::plot_annotation(theme = theme(plot.title = element_text(face = "bold", size = 20)))

ggsave(paste0(fig.path,'QCplot.pdf'),dpi = 300, height = 3, width = 8)

## Batch effects----
p <- DimPlot(preBatch, group.by = 'orig.ident', pt.size = 1) +
  scale_color_discrete(limits = c('Young1','Elderly1','Elderly2')) +
  theme_void() + theme(aspect.ratio = 1, plot.title = element_blank(), text = element_text(size = 12))
p + theme(legend.position = 'none')
ggsave(paste0(fig.path,'preBatch_UMAP.pdf'),dpi = 300, height = 4, width = 4)
ggpubr::as_ggplot(ggpubr::get_legend(p))


p <- DimPlot(postClust, group.by = 'orig.ident', pt.size = 1) +
  scale_color_discrete(limits = c('Young1','Elderly1','Elderly2')) +
  theme_void() + theme(aspect.ratio = 1, plot.title = element_blank(), text = element_text(size = 12))
p + theme(legend.position = 'none')
ggsave(paste0(fig.path,'postBatch_UMAP.pdf'),dpi = 300, height = 4, width = 4)

## dot plot(preAnnot)----
se <- postClust
marker.genes <- list(
  'Steroidogenic' = c('NR5A1','CYP11A1','CYP11B2','CYP11B1','HSD3B2','CYP17A1','CYP21A2'),
  "ZG" = c("DACH1","VSNL1"),
  'ZG/ZF' = c('CCN3','NCAM1'),
  "ZR" = c("CYB5A","SULT2A1"),
  "Cap" = c('NR2F2','RSPO3'),
  "Med" = c('TH','CHGA'),
  "Endo" = c('PECAM1','EMCN'),
  "Fib" = c('LUM','PDGFRA'),
  "VSM" = c("ACTA2","MYH11"),
  "Adipo" = c("ADIPOQ","FABP4"),
  "Mac" = c('CD68','CSF1R'),
  "Tcell" = c('CD4','TRBC1'),
  "Bcell" = c('CD38',"SDC1"),
  'Neut' = c('S100A8','S100A9')
)
DotPlot(se, features = marker.genes, group.by = 'seurat_clusters', cols = c("lightgrey","red"),scale = T,
        col.min = -1, col.max = 1, dot.min = 0.1) + 
  scale_color_gradientn(name = 'Scaled\nexpression', 
                        colors = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)) +
  theme_cowplot() + ylab('Clusters') +
  guides(color = guide_colorbar(title = 'Mean expression\nin group (Z-score)', barwidth = unit(1.5,'in'),
                                direction = 'horizontal', title.position = 'top', title.hjust = 0.5), 
         size = guide_legend(title = 'Fraction of cells\nin group (%)',label.position = 'bottom',keywidth = unit(0.2,'in'),
                             direction = 'horizontal', title.position = 'top', title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1), axis.title.x = element_blank())
ggsave(paste0(fig.path,'preannotDotplot.pdf'), height = 5, width = 14, dpi = 300)


## UMAP annotation----
se <- postAnnot
DimPlot(se,group.by = 'annotation',label=T, pt.size = 1) + 
  scale_color_manual(values = annot_pal) + 
  theme_void() + theme(aspect.ratio = 1, legend.position = 'none', title = element_blank())
ggsave(paste0(fig.path,'UMAP_annot.pdf'), height = 4, width = 4, dpi = 300)

## dot plot(postAnnot)----
se <- postAnnot
marker.genes <- list(
  'Steroidogenic' = c('NR5A1','CYP11A1','CYP11B2','CYP11B1','HSD3B2','CYP17A1','CYP21A2'),
  "ZG" = c("DACH1","VSNL1"),
  'ZG/ZF' = c('CCN3','NCAM1'),
  "ZR" = c("CYB5A","SULT2A1"),
  "Cap" = c('COL1A2','RSPO3'),
  "Med" = c('TH','PNMT'),
  "Endo" = c('PECAM1','EMCN'),
  "Fib" = c('LUM','PDGFRA'),
  "VSM" = c("ACTA2","MYH11"),
  "Adipo" = c("ADIPOQ","FABP4"),
  "Mac" = c('CD68','CSF1R'),
  "Tcell" = c('CD4','TRBC1'),
  "Bcell" = c('CD38',"SDC1"),
  'Neut' = c('S100A8','S100A9')
)

DotPlot(se, features = marker.genes, group.by = 'annotation', cols = c("lightgrey","red"),scale = T,
          col.min = -1, col.max = 1, dot.min = 0.1) + 
    scale_color_gradientn(name = 'Scaled\nexpression', 
                          colors = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)) +
    theme_cowplot() + ylab('Annotation') +
    guides(color = guide_colorbar(title = 'Mean expression\nin group (Z-score)', barwidth = unit(1.5,'in'),
                                  direction = 'horizontal', title.position = 'top', title.hjust = 0.5), 
           size = guide_legend(title = 'Fraction of cells\nin group (%)',label.position = 'bottom',keywidth = unit(0.2,'in'),
                               direction = 'horizontal', title.position = 'top', title.hjust = 0.5)) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1), axis.title.x = element_blank())

ggsave(paste0(fig.path,'postannotDotplot.pdf'), height = 5, width = 14, dpi = 300)

## Bar plot(cell composition)----
se <- postAnnot
se@meta.data %>% as.data.frame() %>% 
  filter(annotation %in% c('ZG','ZF','ZR')) %>% 
  mutate(annotation = factor(annotation,levels = c('ZG','ZF','ZR'))) %>%
  mutate(Age = factor(Age, levels = c('Young', 'Elderly'))) %>%
  group_by(Age,annotation) %>%
  summarise(count = n()) %>% ungroup() %>%
  ggplot(aes(x = Age, y = count)) +
  geom_bar(aes(fill = annotation), color = NA, width = 0.7,alpha = 0.9, 
           position = 'fill', stat = 'identity') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = annot_pal, name = '') +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(), aspect.ratio = 2,
        text = element_text(size = 12),
        axis.title = element_blank(), legend.position = 'right')
ggsave(paste0(fig.path,'cellcomposition.pdf'), height = 3.5, width = 3.5, dpi = 300)


## Heatmap(Elderly vs Young)----
# Heatmap
se <- postAnnot
steroidogenic_gene <- c('CYP11A1','CYP11B2','CYP11B1','HSD3B2','CYP17A1','CYP21A2','CYB5A','STAR')
df <- lapply(c('ZG','ZF','ZR'),function(i){
  sub_se <- subset(se,annotation==i)
  sub_se <- PrepSCTFindMarkers(sub_se)
  Idents(sub_se) <- sub_se$Age
  DEG <- FindMarkers(sub_se[steroidogenic_gene,],assay = 'SCT',
                     # test.use = 'MAST',latent.vars = 'orig.ident',recorrect_umi = F,
                     ident.1 = 'Elderly',ident.2 = 'Young',
                     min.pct = 0,logfc.threshold = 0,
                     densify = T)
  DEG$cluster <- i
  DEG <- DEG %>% rownames_to_column(var = 'gene')
  return(DEG)
}) %>% bind_rows()
## matrix(logFC)
mat <- df %>% dplyr::select(gene,avg_log2FC,cluster) %>% 
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>% 
  column_to_rownames(var = 'gene') %>% as.matrix()
## matrix(pval)
mat2 <- df %>% dplyr::select(gene,p_val_adj,cluster) %>% 
  pivot_wider(names_from = cluster, values_from = p_val_adj) %>% 
  column_to_rownames(var = 'gene') %>% as.matrix()
mat[is.na(mat)] <- 0 ; mat2[is.na(mat2)] <- 0

col_fun <- colorRamp2(breaks = seq(from = -1, to = 1, length = 100), 
                      colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
ht <- Heatmap(mat,
              col = col_fun,column_names_rot = 90, column_names_centered = F,
              heatmap_legend_param = list(title = "log2FC"),cluster_columns = F, 
              width = unit(0.4*ncol(mat),'in'), height = unit(0.4*nrow(mat),'in'),
              # DEGに●を表示
              cell_fun = function(j, i, x, y, w, h, fill){
                grid.text(ifelse(abs(mat[i,j]) > 0.25 & mat2[i,j] < 0.05, "●",''),x,y)}
)
draw(ht)
p <- grid.grabExpr(draw(ht)) 
ggsave(p, filename = paste0(fig.path,'steroidHeatmap.pdf'), height = 5, width = 4, dpi = 300)

## GSEA dot plot----
# plot 
df <- gsea_res %>% 
  mutate(pathway = gsub(pathway,pattern = 'HALLMARK_',replacement = '')) %>%
  mutate(pathway = gsub(pathway,pattern = 'GOBP_',replacement = '')) %>%
  mutate(pathway = gsub(pathway,pattern = '_',replacement = ' ')) %>%
  mutate(pathway = str_to_title(pathway)) 

hclust_df <- df %>%
  # filter(padj < 0.05) %>%
  dplyr::select(cluster,pathway,NES) %>% 
  pivot_wider(values_from = NES, names_from = cluster) %>% column_to_rownames(var = 'pathway')
clust <- hclust(dist(hclust_df %>% as.matrix())) # hclust with distance matrix

dotplot <- df %>% 
  mutate(cluster = factor(cluster, levels = levels(se$annotation))) %>% 
  mutate(pathway = factor(pathway,levels = clust$labels[clust$order])) %>% # hclust order
  filter(padj < 0.05) %>%
  mutate(NES = ifelse(NES > 2, 2, NES)) %>% 
  mutate(NES = ifelse(NES < -2, -2, NES)) %>% 
  ggplot(aes(x = cluster, y = pathway)) +
  geom_point(aes(color = NES, size = -log10(padj)),alpha = 1) +
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100), limits = c(-2,2)) +
  theme_bw() +
  theme(axis.line  = element_blank(), aspect.ratio = 3,
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(position = "right") +
  theme(axis.ticks = element_blank())
dotplot
ggsave(filename = paste0(fig.path,'GSEA_dot_hallmark.pdf'), height = 6, width = 5.5, dpi = 300)

## RSS plot----
for(x in unique(rss$annotation)){
  tmp <- rss %>% filter(annotation == x)
  label_df <- tmp %>% filter(rank <= 5)
  tmp %>% 
    ggplot(aes(x = rank, y = RSS)) +
    geom_point(color = "lightgray", size = 2, alpha = 0.5) +
    geom_point(data = label_df,aes(x = rank, y = RSS), 
               color = "black", size = 2, alpha = 1) +
    ggrepel::geom_label_repel(data = label_df,aes(x = rank, y = RSS, label = TF), size = 2.5,
                              box.padding = 0.5, max.overlaps = Inf) +
    ggtitle(x) +
    ylim(c(NA,max(tmp$RSS)*1.1)) +
    theme_classic() +
    theme(aspect.ratio = 1.5)
  ggsave(paste0(fig.path,'RSSplot_',x,'.pdf'), height = 3, width = 2.2, dpi = 300)
}

## Venn diagram(TFs common to young and elderly)----
TF_list <- lapply(unique(rss$annotation),function(x){
  rss %>% filter(annotation == x) %>% top_n(n = 10, wt = -rank) %>% .$TF
})
names(TF_list) <- unique(rss$annotation)
for(annotation in c('ZG','ZF','ZR')){
  gene_df <- data.frame(gene = unique(c(TF_list[[paste0('Young_',annotation)]],TF_list[[paste0('Elderly_',annotation)]])),
                        row.names = 'gene',
                        set1 = unique(c(TF_list[[paste0('Young_',annotation)]],TF_list[[paste0('Elderly_',annotation)]])) %in% TF_list[[paste0('Young_',annotation)]],
                        set2 = unique(c(TF_list[[paste0('Young_',annotation)]],TF_list[[paste0('Elderly_',annotation)]])) %in% TF_list[[paste0('Elderly_',annotation)]])
  colnames(gene_df) <- c('Young','Elderly')
  
  ggvenn::ggvenn(gene_df, text_size = 5, set_name_size = F,
                 columns = colnames(gene_df),
                 show_percentage = F, fill_alpha = .4,
                 # auto_scale = T,
                 stroke_color = 'grey',
                 fill_color = unname(age_pal[colnames(gene_df)])) +
    theme(aspect.ratio = 0.65, text = element_text(size = 12))
  ggsave(paste0(fig.path,'Venn_TFtop10_',annotation,'.pdf'), dpi = 300, height = 1.5, width = 2)  
}


## Trajectory plot----
se <- Traj

reduction <- 'dc' ; dims <- c(1,2)
DimPlot(se, reduction = reduction, dims = dims, group.by = 'annotation', label = F, pt.size = 1.5) +
  scale_color_manual(values = annot_pal) +
  theme_void() + theme(aspect.ratio = 1, legend.position = 'none', title = element_blank())
ggsave(paste0(fig.path,'Diffusion_annot.pdf'),dpi = 300, height = 4, width = 4)

# trajectory
p <- FeaturePlot(se, reduction = reduction, dims = dims, features = 'Pseudotime', label = F, pt.size = 1.5) +
  scale_color_viridis_c(name = 'Pseudotime')
##legend
ggpubr::as_ggplot(ggpubr::get_legend(p))
ggsave(filename = paste0(fig.path,'Pseudotime_legend.pdf'), height = 4, width = 1.5, dpi = 300)

p + geom_path(data = data.frame(slingCurves(pto)[[1]]$s[slingCurves(pto)[[1]]$ord, ]) %>% dplyr::select(xaxis = dims[1], yaxis = dims[2]),
              aes(x = xaxis, y = yaxis), colour = "grey", linewidth = 1, linetype = 1, arrow = grid::arrow(ends = 'last', length = unit(x = 5,units = 'mm')))+
  theme_void() + theme(aspect.ratio = 1, legend.position = 'none', title = element_blank())
ggsave(paste0(fig.path,'Diffusion_pseudotime.pdf'),dpi = 300, height = 4, width = 4)


## WNT/PKA along pseudotime----
se <- ssGSEA_Young
data.frame(
  WNT = as.numeric(se@assays$ssgsea['HALLMARK-WNT-BETA-CATENIN-SIGNALING',]),
  PKA = as.numeric(se@assays$ssgsea['GOBP-PROTEIN-KINASE-A-SIGNALING',]),
  Annot = se$annotation,
  Pseudotime = se$Pseudotime,
  Age = factor(se$Age,levels = c('Young','Elderly'))
  ) %>% 
  dplyr::select(Pseudotime, Age, ES = PKA) %>% 
  rstatix::cor_test(Pseudotime,ES)

ggplot(aes(x = Pseudotime, y = ES)) +
  geom_point(aes(color = Annot), alpha = 0.8) +
  geom_smooth(method = 'lm',se = F) +
  scale_color_manual(values = annot_pal) +
  labs(y = 'Enrichment score', x = 'Pseudotime') +
  theme_bw() + 
  theme(aspect.ratio = 1, legend.title = element_blank(),legend.position = 'none', 
        axis.title = element_text(color = 'black', size = 12), 
        axis.text = element_text(color = 'black', size = 12),
        text = element_text(color = 'black', size = 12))

## LR alluvial ZF and Macrophage----
eat.pathways <- readRDS("eat_pathways.rds") 
cols = c('Find_me_signalling' = '#E41A1C',	'Eat_me_signalling' = '#377EB8',
         'Dont_eat_me' = '#4DAF4A','Phagocytosis_signalling' = '#984EA3', 'other' = 'black')

lrobj_tbl <- LRobj@tables$EXP_x_CTR %>% 
  filter(Ligand %in% unlist(eat.pathways) | Receptor %in% unlist(eat.pathways)) %>% # phagocytosisに関連した分子に絞る(DOI: 10.1126/sciadv.add0422)
  arrange(desc(abs(LRScore))) %>% 
  filter(abs(LRScore) > 0.7)  # 有意とした条件

layer <- 'ZF'  
lrobj_tbl %>% 
  filter(Ligand.Cluster == layer & Receptor.Cluster == 'Macrophage') %>% 
  mutate(Freq = 1) %>% 
  mutate(UpDown = factor(ifelse(LRScore>0,'Up','Down'),levels = c('Up','Down'))) %>% 
  dplyr::select(source = Ligand.Cluster, Ligand, Receptor, target = Receptor.Cluster, everything()) %>% 
  to_lodes_form(key = "Demographic",axes = 1:4) %>% 
  mutate(stratum = fct_reorder(stratum, -LRScore)) %>%
  mutate(Demographic = factor(Demographic, levels = c("source", "Ligand", "Receptor", "target"))) %>% 
  
  mutate(Signal = ifelse(stratum %in% eat.pathways$Find_me_signalling, "Find_me_signalling",
                         ifelse(stratum %in% eat.pathways$Eat_me_signalling, "Eat_me_signalling",
                                ifelse(stratum %in% eat.pathways$Dont_eat_me, "Dont_eat_me",  
                                       ifelse(stratum %in% eat.pathways$Phagocytosis_signalling, "Phagocytosis_signalling","other"))))) %>% 
  
  ggplot(aes(x = Demographic, stratum = stratum, alluvium = alluvium, y = Freq, label = stratum)) +
  ggalluvial::geom_alluvium(aes(fill = UpDown),width = 1/12, discern = FALSE) +
  ggalluvial::geom_stratum(width = 1/12) + 
  geom_label(stat = 'stratum',size = 4, aes(color = Signal)) +
  scale_x_discrete(expand = c(0.05, 0.05),) +
  scale_y_continuous(breaks = NULL) + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), 
        aspect.ratio = 1, panel.grid.major = element_blank(),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(name = 'LR Score', values = c('#D73027','#4575B4'))
  
ggsave(filename = paste0(fig.path,'LRanalysis_',layer,'toMacrophage','.pdf'),dpi = 300, height = 5, width = 7)
  

## Macrophage GSEA bar plot----
select_pathway1 <- c('HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_BILE_ACID_METABOLISM',
                     'HALLMARK_ADIPOGENESIS','HALLMARK_FATTY_ACID_METABOLISM','M2_macrophage')
select_pathway2 <- c('HALLMARK_ALLOGRAFT_REJECTION','HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                     'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_KRAS_SIGNALING_UP',
                     'M1_macrophage')
select_pathway <- c(select_pathway1,select_pathway2)

GSEA_Mac %>% 
  mutate(Age = ifelse(NES > 0,'Elderly','Young')) %>% 
  filter(pathway %in% select_pathway) %>% 
  mutate(pathway = gsub(pathway, pattern = 'HALLMARK_', replacement = '')) %>%
  mutate(pathway = gsub(pathway,pattern = '_', replacement = ' ')) %>% 
  mutate(pathway = str_to_title(pathway)) %>% 
  ggplot(aes(x = NES, y = fct_reorder(pathway,NES), fill = Age)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = age_pal) +
  scale_y_discrete(labels = scales::wrap_format(80)) +
  xlim(c(-max(abs(gsea_res$NES)),max(abs(gsea_res$NES)))) +
  theme_bw() + 
  theme(
    aspect.ratio = 1, 
    legend.position = 'none',
    panel.grid.minor.y =  element_blank(), 
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12,color = 'black'),
    axis.text.y = element_text(size = 12,color = 'black'),
    axis.title.x = element_text(size = 15,color = 'black'),
    legend.text = element_text(size = 12,color = 'black'), 
    legend.title = element_text(size = 12,color = 'black')
  )

ggsave(paste0(fig.path,'GSEA_Mac_barplot.pdf'),height = 4, width = 7, dpi = 300)

# ST------------------------------------------------------------
## setup----
# load library
easypackages::libraries(c('Seurat','tidyverse','patchwork','harmony',
                          'ComplexHeatmap','circlize','RColorBrewer', 'cowplot'))

# load objects
save.path <- 'object folder path'
preQC <- readRDS(paste0(save.path,'preQC.rds'))
postQC <- readRDS(paste0(save.path,'postQC.rds'))
preBatch<- readRDS(paste0(save.path,'preBatch.rds'))
postClust <- readRDS(paste0(save.path,'postClust.rds'))
postAnnot <- readRDS(paste0(save.path,'postAnnot.rds'))
assort_res <- readRDS(paste0(save.path,'assort_res.rds'))
Deconv <- readRDS(paste0(save.path,"Deconv.rds"))
ssGSEA <- readRDS(paste0(save.path,"ssGSEA.rds"))

# color palette
age_pal <- c('#C77CFF','#7CAE00') ; names(age_pal) <- c('Elderly','Young')
sex_pal <- c('#00BFC4','#F8766D') ; names(sex_pal) <- c('Male','Female')
annot_pal_SRT <- c('Capsule' = '#66C2A5', 'ZG' = '#FC8D62', 'ZF' = '#8DA0CB', 'ZR' = '#E78AC3',
                   'Medulla' = '#A6D854', 'Stromal' = '#FFD92F')


# figure save folder
fig.path <- 'figure folder path'

## QCplot----
plot.df <- lapply(preQC,function(se){return(se@meta.data)}) %>% bind_rows()
plot.df$orig.ident <- factor(plot.df$orig.ident, levels = c('Young1','Elderly1','Elderly2'))
plot.set <- theme_classic() + 
  theme(aspect.ratio = 1, axis.title = element_blank(), 
        axis.text.x = element_text(size = 10), axis.ticks.x = element_blank(), 
        title = element_text(face = "bold"), legend.position = "right")
p1 <- plot.df %>% ggplot(aes(x = orig.ident, y = nCount_Spatial)) +
  labs(title = "#UMIs") + geom_violin(color = "black", fill = "lightgrey") +
  plot.set
p2 <- plot.df %>% ggplot(aes(x = orig.ident, y = nFeature_Spatial)) + 
  labs(title = "#Genes") + 
  geom_violin(color = "black", fill = "lightgrey") + 
  geom_hline(yintercept = feature.min, linetype = "dashed", color = "red") +
  plot.set
wrap_plots(p1,p2) + 
  patchwork::plot_annotation(theme = theme(plot.title = element_text(face = "bold", size = 20)))
ggsave(filename = paste0(fig.path,'QCplot.pdf'),dpi = 300, height = 3, width = 6)


# Spatial plot
names(preQC) <- c('Young1','Elderly1','Elderly2')
for(sample in c('Young1','Elderly1','Elderly2')){
  # HE
  se <- se.list[[sample]]
  p <- SpatialDimPlot(se, images = sample, alpha = 0, crop = F) & theme(legend.position = "none")
  ggsave(plot = p, filename = paste0(fig.path,'QCplot_HE_',sample,'.pdf'), dpi = 300, height = 4, width = 4)
  
  # UMI
  breaks <- floor(max(se$nCount_Spatial)/2/10000)*10000
  p <- SpatialFeaturePlot(se,images = sample, crop = F,features = c("nCount_Spatial"), pt.size.factor = 1.3, stroke = NA) & 
    scale_fill_viridis_c(name = '#UMIs', option = 'magma', breaks = c(0,breaks,breaks*2)) &
    theme(legend.position = "right")
  ggsave(plot = p + theme(legend.position = 'none'), 
         filename = paste0(fig.path,'QCplot_UMI_',sample,'.pdf'), dpi = 300, height = 4, width = 4)

  # gene count
  breaks <- floor(max(se$nFeature_Spatial)/2/1000)*1000
  p <- SpatialFeaturePlot(se, images = sample, crop = F,features = c("nFeature_Spatial"), pt.size.factor = 1.3, stroke = NA) & 
    scale_fill_viridis_c(name = '#Genes', option = 'viridis', breaks = c(0,breaks,breaks*2)) &
    theme(legend.position = "right")
  ggsave(plot = p + theme(legend.position = 'none'), 
         filename = paste0(fig.path,'QCplot_Gene_',sample,'.pdf'), dpi = 300, height = 4, width = 4)

  # Exclude
  se@meta.data <- se@meta.data %>%
    mutate(QC = ifelse(nFeature_Spatial>feature.min & Histological_annotation != "Exclude",NA,'Exclude'))
  p <- SpatialDimPlot(se,images = sample, crop = F,group.by = c("QC"),alpha = 1, pt.size.factor = 1.3,stroke = NA) &
    theme(legend.position = "none")
  ggsave(plot = p, filename = paste0(fig.path,'QCplot_Exclude_',sample,'.pdf'), dpi = 300, height = 4, width = 4)
  
}

## Batch effects----
p <- DimPlot(preBatch, group.by = 'orig.ident', pt.size = 1) +
  scale_color_discrete(limits = c('Young1','Elderly1','Elderly2')) +
  theme_void() + theme(aspect.ratio = 1, plot.title = element_blank(), text = element_text(size = 12))
p + theme(legend.position = 'none')
ggsave(paste0(fig.path,'preBatch_UMAP.pdf'),dpi = 300, height = 4, width = 4)
ggpubr::as_ggplot(ggpubr::get_legend(p))


p <- DimPlot(postClust, group.by = 'orig.ident', pt.size = 1) +
  scale_color_discrete(limits = c('Young1','Elderly1','Elderly2')) +
  theme_void() + theme(aspect.ratio = 1, plot.title = element_blank(), text = element_text(size = 12))
p + theme(legend.position = 'none')
ggsave(paste0(fig.path,'postBatch_UMAP.pdf'),dpi = 300, height = 4, width = 4)

## UMAP annotation(preAnnot)----
se <- postClust
DimPlot(se,group.by = 'seurat_clusters',label=T, pt.size = 1) + 
  theme_void() + theme(aspect.ratio = 1, legend.position = 'none', title = element_blank())
ggsave(paste0(fig.path,'UMAP_preannot.pdf'), height = 4, width = 4, dpi = 300)

## spatial plot(preAnnot)----
se <- postClust
for(i in unique(se$orig.ident)){
  SpatialDimPlot(se, images = i, group.by = 'seurat_clusters', crop = F, stroke = NA, pt.size.factor = 1.3) +
    theme_void() + theme(legend.position = 'none')
  ggsave(filename = paste0(fig.path,"preAnnotSpatial_",sample,'.pdf'), height = 4, width = 4, dpi = 300)
}

## dot plot(preAnnot)----
se <- postClust
marker.genes <- list(
  'Steroidogenic' = c('NR5A1','CYP11A1','CYP11B2','CYP11B1','CYP17A1','CYP21A2'),
  "ZG" = c("DACH1","VSNL1"),
  'ZG/ZF' = c('CCN3','NCAM1'),
  "ZR" = c("CYB5A","SULT2A1"),
  "Cap" = c('NR2F2','RSPO3'),
  "Med" = c('TH','CHGA'),
  "Endo" = c('PECAM1','EMCN'),
  "Fib" = c('LUM','PDGFRA'),
  "VSM" = c("ACTA2","MYH11"),
  "Adipo" = c("ADIPOQ","FABP4"),
  "Mac" = c('CD68','CSF1R'),
  "Tcell" = c('CD4','TRBC1'),
  "Bcell" = c('CD38',"SDC1"),
  'Neut' = c('S100A8','S100A9')
)
DotPlot(se, features = marker.genes, group.by = 'seurat_clusters', cols = c("lightgrey","red"),scale = T,
        col.min = -1, col.max = 1, dot.min = 0.1) + 
  scale_color_gradientn(name = 'Scaled\nexpression', 
                        colors = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)) +
  theme_cowplot() + ylab('Clusters') +
  guides(color = guide_colorbar(title = 'Mean expression\nin group (Z-score)', barwidth = unit(1.5,'in'),
                                direction = 'horizontal', title.position = 'top', title.hjust = 0.5), 
         size = guide_legend(title = 'Fraction of spots\nin group (%)',label.position = 'bottom',keywidth = unit(0.2,'in'),
                             direction = 'horizontal', title.position = 'top', title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1), axis.title.x = element_blank())
ggsave(paste0(fig.path,'preannotDotplot.pdf'), height = 5, width = 14, dpi = 300)

## UMAP annotation----
se <- postAnnot
DimPlot(se,group.by = 'annotation',label=T, pt.size = 1) + 
  scale_color_manual(values = annot_pal) + 
  theme_void() + theme(aspect.ratio = 1, legend.position = 'none', title = element_blank())
ggsave(paste0(fig.path,'UMAP_annot.pdf'), height = 4, width = 4, dpi = 300)

## spatial plot(postAnnot)----
se <- postClust
for(i in unique(se$orig.ident)){
  SpatialDimPlot(se, images = i, group.by = 'annotation', cols = annot_pal_SRT, crop = F, stroke = NA, pt.size.factor = 1.3) +
    theme_void() + theme(legend.position = 'none')
  ggsave(filename = paste0(fig.path,"postAnnotSpatial_",sample,'.pdf'), height = 4, width = 4, dpi = 300)
}

## dot plot(postAnnot)----
se <- postAnnot

DotPlot(se, features = marker.genes, group.by = 'annotation', cols = c("lightgrey","red"),scale = T,
        col.min = -1, col.max = 1, dot.min = 0.1) + 
  scale_color_gradientn(name = 'Scaled\nexpression', 
                        colors = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)) +
  theme_cowplot() + ylab('Annotation') +
  guides(color = guide_colorbar(title = 'Mean expression\nin group (Z-score)', barwidth = unit(1.5,'in'),
                                direction = 'horizontal', title.position = 'top', title.hjust = 0.5), 
         size = guide_legend(title = 'Fraction of spots\nin group (%)',label.position = 'bottom',keywidth = unit(0.2,'in'),
                             direction = 'horizontal', title.position = 'top', title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1), axis.title.x = element_blank())

ggsave(paste0(fig.path,'postannotDotplot.pdf'), height = 5, width = 14, dpi = 300)

## Bar plot(cell composition)----
se <- postAnnot
se@meta.data %>% as.data.frame() %>% 
  filter(annotation %in% c('ZG','ZF','ZR')) %>% 
  mutate(annotation = factor(annotation,levels = c('ZG','ZF','ZR'))) %>%
  mutate(Age = factor(Age, levels = c('Young', 'Elderly'))) %>%
  group_by(Age,annotation) %>%
  summarise(count = n()) %>% ungroup() %>%
  ggplot(aes(x = Age, y = count)) +
  geom_bar(aes(fill = annotation), color = NA, width = 0.7,alpha = 0.9, 
           position = 'fill', stat = 'identity') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = annot_pal, name = '') +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(), aspect.ratio = 2,
        text = element_text(size = 12),
        axis.title = element_blank(), legend.position = 'right')
ggsave(paste0(fig.path,'cellcomposition.pdf'), height = 3.5, width = 3.5, dpi = 300)

## highlight layers----
se <- postAnnot
for(sample in unique(se$orig.ident)){
  for(layer in c('ZG','ZF','ZR')){
    SpatialDimPlot(se, alpha = NULL, images = sample, stroke = NA, crop = F, pt.size.factor = 1.3,
                        cols.highlight = c(annot_pal_SRT[[layer]], alpha("#f5f5f5", 0)),
                        cells.highlight = Cells(se)[se$annotation==layer]) +
      theme_void() + theme(plot.title = element_blank(), legend.position = 'none')
    ggsave(filename = paste0(fig.path,'3layer_',sample,'_',layer,'.pdf'),width = 4, height = 4, dpi = 300)
  }}

## ROIs------------
spata.list <- lapply(c('Young1','Elderly1','Elderly2'),function(sample){
  loadSpataObject(directory_spata = paste0(save.path,'SPATA_',sample,'.rds'))
})
names(spata.list) <- c('Young1','Elderly1','Elderly2')

for(sample in c('Young1','Elderly1','Elderly2')){
  for(id in c('ROI1','ROI2','ROI3','ROI4')){
    ROI <- getTrajectory(spata.list[[sample]],id = id)@projection %>% pull(barcodes)
    df <- data.frame(barcodes = getBarcodes(spata.list[[sample]]),group = 'out') %>%
      mutate(group = ifelse(barcodes %in% ROI, 'in', group))
    tmp <- addFeatures(object = spata.list[[sample]], feature_df = df, overwrite = TRUE)
    outline <- ggpLayerGroupOutline(object = tmp, grouping = "group", groups_subset = "in", line_size = 0.8) 
    
    plotSpatialTrajectories(spata.list[[sample]], ids = id, color_by = 'annotation', 
                            display_image = T, sgmt_size = NA, pt_alpha = 0.5, pt_alpha2 = 1, pt_size = 1.5) +
      scale_color_manual(values = annot_pal_SRT) + theme(legend.position = 'none') + outline
    ggsave(filename = paste0(fig.path,'ROIs_',sample,'.pdf'), height = 4, width = 16, dpi = 300)
  }
}

## Assortativity box plot----
assort_res %>% 
  filter(label %in% c('ZG','ZF','ZR')) %>% 
  mutate(Age = factor(ifelse(Sample=='Young1','Young','Elderly'),levels = c('Young','Elderly'))) %>% 
  group_by(label,Age) %>% 
  summarize(
    mean = mean(avg_k_scaled,na.rm=TRUE),
    sd = sd(avg_k_scaled,na.rm=TRUE)
  ) %>% 
  ggplot(aes(x = Age, y = mean)) +
  geom_boxplot(aes(fill = Age), alpha = 0.6, outlier.alpha = 1, outlier.size = .3,
               fatten = 1.5, size = 0.3, width = 0.5) + 
  scale_fill_manual(values = age_pal) +
  ylim(c(NA,0.9)) +
  ylab("Assortativity score") +
  facet_wrap(~label) + cowplot::theme_cowplot() +
  theme(aspect.ratio = 3,
        plot.title = element_blank(), 
        legend.position = 'none',
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 15))

ggsave(paste0(fig.path,'Assortativity.pdf'),height = 5, width = 5, dpi = 300)

## Deconvolution Macrophage spatial plot----
se <- Deconv
for(i in unique(se$orig.ident)){
  SpatialFeaturePlot(se %>% subset(annotation %in% c('ZG','ZF','ZR')),
                     images = i,features = 'Prob_Macrophage',stroke = NA, crop = F, pt.size.factor = 1.3) & 
    scale_fill_gradient(low = alpha('lightgrey',alpha = 0.1), high = 'red', name = NULL,
                        limits = c(0.0,0.3), breaks = c(0.0, 0.15, 0.3)) &
    theme(title = element_blank(), legend.position = 'none')
  ggsave(paste0(fig.path,'Deconv_Macrophage_',i,'.pdf'), dpi = 300, height = 4, width = 4)
}

## Deconvolution Macrophage box plot----
se <- Deconv

# Statistical tests for macrophage proportions
## Young vs Elderly
young <- se$Prob_Macrophage[se$Age == 'Young']
elderly <- se$Prob_Macrophage[se$Age == 'Elderly']
lapply(list(young,elderly),summary)
wilcox.test(young,elderly)
## Elderly each layer
zg <- se$Prob_Macrophage[se$Age == 'Elderly' & se$annotation == 'ZG']
zf <- se$Prob_Macrophage[se$Age == 'Elderly' & se$annotation == 'ZF']
zr <- se$Prob_Macrophage[se$Age == 'Elderly' & se$annotation == 'ZR']
kruskal.test(list(zg,zf,zr))
lapply(list(zg,zf,zr),summary)
wilcox.test(zg,zf)$p.value %>% p.adjust(.,method = 'BH', n = 3)
wilcox.test(zg,zr)$p.value %>% p.adjust(.,method = 'BH', n = 3)
wilcox.test(zf,zr)$p.value %>% p.adjust(.,method = 'BH', n = 3)

tibble(Proportion = se$Prob_Macrophage, Layer = se$annotation, Sample = se$orig.ident) %>% 
  filter(Layer %in% c('ZG','ZF','ZR')) %>% 
  mutate(Sample = factor(ifelse(Sample=='Young1','Young','Elderly'),levels = c('Young','Elderly'))) %>% 
  ggplot(aes(x = Layer, y = Proportion, fill = Layer)) +
  geom_boxplot(aes(fill = Age), 
               outlier.alpha = 1, outlier.size = .5,
               fatten = 1.5, size = 0.5, width = 0.5,
               alpha = 0.6) + 
  scale_fill_manual(values = annot_pal) +
  facet_wrap(~Sample) + cowplot::theme_cowplot(font_size = 12) + 
  theme(aspect.ratio = 1.5, text = element_text(size = 12), 
        axis.title.x = element_blank(), legend.position = 'none')
ggsave(paste0(fig.path,'Deconvolution_violin_Macrophage.pdf'),dpi = 300, height = 4, width = 5.5)

## correlation Macrophage proportion and SenMayo----
se <- ssGSEA

cor_df <- apply(GetAssayData(se,assay = 'ssgsea',slot = 'data'),MARGIN = 1,
                FUN = function(x){
                  tmp <- cor.test(x,se$Prob_Macrophage)
                  r <- tmp$estimate ; p <- tmp$p.value
                  data.frame(r = r, p = p)
                }) %>% bind_rows()
cor_df$pathway <- rownames(se@assays$ssgsea)
cor_df %>% filter(p < 0.05) %>% arrange(-r)

data.frame(SenMayo = as.numeric(se@assays$ssgsea['SenMayo',]),
           INF = as.numeric(se@assays$ssgsea['HALLMARK-INFLAMMATORY-RESPONSE',]),
           Macrophage = se$Prob_Macrophage,
           Annot = se$annotation,
           Age = se$Age) %>% 
  filter(Annot %in% c('ZG','ZF','ZR')) %>%
  # rstatix::cor_test(SenMayo,Macrophage,method = 'pearson')
  ggplot(aes(x = SenMayo, y = Macrophage)) +
  geom_point(aes(color = Annot)) + geom_smooth(method = 'lm', se = F) +
  scale_color_manual(values = annot_pal_SRT) + 
  labs(x = 'NES', y = 'Macrophage proportion') +
  theme_bw() + theme(aspect.ratio = 1, legend.title = element_blank(), text = element_text(size = 12))

ggsave(filename = paste0(fig.path,'Macrophage_SenMayo.pdf'), width = 4, height = 3, dpi = 300)


