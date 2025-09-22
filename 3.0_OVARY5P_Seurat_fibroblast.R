library(ggrepel)

####Fibroblast####
wd.path <- "d:/project_directory_name/ovary_project/Seurat/final/"
fibro.path <- paste0(wd.path, "stromal/fibroblast/")
setwd(fibro.path)

fibro <- subset(ovary, idents = "Fibroblast")

table(fibro$orig.ident, fibro$chemo)


Idents(fibro) <- "chemo"
clustering <- fibro@active.ident
cellsIn <- names(clustering[clustering == "Post-treatment"])
cellsOut <- names(clustering[clustering == "Pre-treatment"])
fibro<- ScaleData(fibro)

fibro_markers <- FindMarkers(fibro, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)

fibro_markers$significant <- ifelse(-log10(fibro_markers$p_val_adj) > 5 & abs(fibro_markers$avg_log2FC) >= 1, 
                                    ifelse(fibro_markers$avg_log2FC >= 1 , 'Up' ,'Down'),' ')

fibro_markers$significant <- fibro_markers$significant %>% as.factor()
fibro_markers$significant %>% table

fibro_markers$significant <- factor(fibro_markers$significant, levels = c("Up", "Down", " "))

ggplot(fibro_markers, aes(avg_log2FC, -log10(p_val_adj), fill = significant, color = significant)) + 
  geom_point(size = 1.5, shape = 21) +
  theme_bw() + 
  geom_hline(yintercept=5, linetype='dashed', color='red', size=0.5) +
  geom_vline(xintercept=1, linetype = 'dashed', color='red', size = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color='red', size = 0.5) + ylim(c(0,200)) +
  scale_fill_manual(values=c("red","blue","gray67")) +
  scale_color_manual(values=c("red","blue","gray67")) +
  geom_text_repel(data = fibro_markers[fibro_markers[,"significant"] == "Up" | 
                                         fibro_markers[,"significant"] == "Down" ,],
                  aes(label = genes), 
                  size = 5,
                  color = "black") + theme(legend.text = element_text(size = 15),
                                           legend.title = element_text(size = 15),
                                           axis.text = element_text(size = 15),
                                           axis.title = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA)))) 

ggsave("fibro_markers_volcano.png", width = 8, height = 6, dpi = 300)
ggsave("fibro_markers_volcano.pdf", width = 8, height = 6, dpi = 600)

fibro <- NormalizeData(object=fibro, normalization.method="LogNormalize", scale.factor=10000)

## Identification of highly variable features
fibro <- FindVariableFeatures(fibro, selection.method = "vst", nfeatures = num.features)
top10 <- head(VariableFeatures(fibro), 10)

## Centering the data (z-scoring)
all.genes <- rownames(fibro)
fibro <- ScaleData(fibro, features = all.genes)

## Run PCA
library(tibble)
fibro <- RunPCA(fibro, npcs=50, features = VariableFeatures(fibro))

## Determine statistically significant principal components
#### Elbow plot
ElbowPlot(fibro, ndims=50, reduction="pca") + theme_light()

fibro@reductions$pca@stdev[1:33] %>%  sum() /
  fibro@reductions$pca@stdev %>% sum()

##################################
library(Seurat)
library(bbknnR)
#20, 0.5
fibro <- fibro %>% 
  RunBBKNN(run_TSNE = FALSE, n_pcs = 14, batch_key = "orig.ident") %>% 
  FindClusters(resolution = 0.3, graph.name= "bbknn") 

DimPlot(fibro, group.by = "seurat_clusters")

Idents(fibro) <- "celltype"
fibro_Markers = FindAllMarkers(object=fibro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

fibro_Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
fibro <- ScaleData(fibro, features = top10$gene)

DoHeatmap(object=fibro, features= top10$gene, draw.lines = FALSE, size = 4) + theme(text = element_text(size = 12)) + NoLegend()
ggsave(file="fibro_fb_marker_heatmap.png", width=6, height=8, dpi = 300) 

av_fibro <-AverageExpression(fibro, return.seurat = TRUE)

DoHeatmap(object=av_fibro, features=top10$gene, draw.lines = FALSE, angle = 45, group.colors = celltype.colv,size = 4) + theme(text = element_text(size = 12)) + NoLegend()

ggsave(file="fibro_heatmap.png", width=6.5, height=8, dpi = 300) 
ggsave(file="fibro_heatmap.pdf", width=6.5, height=10, dpi = 600) 


celltype.id <- c("0"="FB0_MMP11",  "1"="FB1_CFD", "2"= "FB2_MYH11", 
                 "3"="FB3_STAR", "4"="FB4_CALB2", "5"="FB5_Doublet")

identical(length(celltype.id), length(unique(fibro$seurat_clusters)))
fibro@meta.data$celltype <- plyr::mapvalues(fibro@meta.data$seurat_clusters, from=names(celltype.id), to=celltype.id)
fibro$celltype %>% table
Idents(fibro) <- "celltype"

celltype.colv <- brewer.pal(n = 8, "Dark2")

FeaturePlot(fibro, c("PDGFC", "THBS1", "COMP", "SKIL"))
FeaturePlot(fibro, c("COMP", "FN1"))
pt.size = 3

p1 <- DimPlot(fibro, reduction = "umap", group.by = "celltype", pt.size= pt.size, cols = celltype.colv) + 
  labs(title = " ")  +
  theme(aspect.ratio= 1, 
        line = element_blank(), 
        rect = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7), ncol = 4, byrow = FALSE))
AugmentPlot(p1)

ggsave("dimplot_celltype_fibro.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("dimplot_celltype_fibro.pdf", width = 5.5, height = 5.5, dpi = 600)

legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_fibro_legend.png", width = 2, height = 2, dpi = 300)
ggsave("dimplot_celltype_fibro_legend.pdf", width = 2, height = 2, dpi = 600)



p2 <- DimPlot(fibro, reduction = "umap", group.by = "orig.ident", pt.size= pt.size, cols = orig.ident.colv) + 
  labs(title = " ") +
  theme(aspect.ratio= 1,
        line = element_blank(), 
        rect = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.title = element_blank())  

AugmentPlot(p2)
ggsave("fibro_dimplot_celltype_orig.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("fibro_dimplot_celltype_orig.pdf", width = 5.5, height = 5.5, dpi = 600)

legend <- get_legend(p2 + theme(legend.background = element_blank(),
                                legend.title = element_blank()) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend, xlim = c(1,1) )
ggsave("fibro_dimplot_celltype_orig_legend.png", width = 1.5, height = 3.5, dpi = 300)
ggsave("fibro_dimplot_celltype_orig_legend.pdf", width = 1.5, height = 3.5, dpi = 600)

fibro$chemo <- factor(fibro$chemo, levels = c("Pre-treatment","Post-treatment"))

p6 <- DimPlot(fibro, reduction = "umap", group.by = "chemo", pt.size=pt.size, cols = c("#00B0F0", "#00B050")) + 
  labs(x = "bbknn_umap_1", y = "bbknn_umap_1", title = " ") +
  theme(aspect.ratio= 1,
        line = element_blank(), 
        rect = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.title = element_blank())
AugmentPlot(p6)
ggsave("fibro_dimplot_total_treatment.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("fibro_dimplot_total_treatment.pdf", width = 5.5, height = 5.5, dpi = 600)



legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("fibro+dimplot_total_treatment_legend.png", width = 3, height = 1, dpi = 300)
ggsave("fibro+dimplot_total_treatment_legend.pdf", width = 3, height = 1, dpi = 600)


Input.df <- fibro@meta.data
Input.df$celltype <- gdata::drop.levels(Input.df$celltype)

Input.df$orig.ident %>% table
Input.df$orig.ident <- factor(Input.df$orig.ident, levels = c("T107", "T70",
                                                              "T04", "T110",
                                                              "T92", "T66",
                                                              "T17", "T18",
                                                              "T39", "T34"))
nCells2 <- dplyr::count (Input.df, orig.ident)

pre_Counted2 <- Input.df %>% group_by(orig.ident) %>% dplyr::count(celltype)
Counted2 <- pre_Counted2 %>% ungroup %>% tidyr::complete(celltype, orig.ident,  fill = list(n = 0)) # with zero-values


Input2 <- merge (x=Counted2, y=nCells2, by = "orig.ident", all.x=TRUE)
Input2$Percent <- Input2$n.x / Input2$n.y * 100

table(Input2$celltype)

count_input <- unique(Input2[,c(1,4)])
Input2$orig.ident %>% table

p1 <- ggplot(Input2, aes(x=orig.ident, y=Percent, fill=celltype)) + 
  geom_bar(stat="identity",position = position_fill(reverse = TRUE), color="black")  +
  labs(x=NULL, y="Proportion") + coord_flip() + scale_fill_manual(values = celltype.colv) +
  theme_classic() + theme(legend.position = "bottom",
                          legend.direction = "horizontal",
                          legend.key.width = unit(0.4,"cm"),
                          legend.key.size =  unit(0.4,"cm"),
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(colour="black", size= 12),
                          axis.text = element_text(size = 12)) + guides(fill=guide_legend(ncol=3))


p1 + ggplot(count_input, aes(x=orig.ident, y=n.y)) + geom_bar(stat="identity")  +
  labs(x=NULL, y="count") + coord_flip() + 
  theme_classic() + theme(aspect.ratio= 2,
                          axis.ticks = element_blank(),
                          legend.position = "right", 
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text = element_text(size = 12))

ggsave("celltype_proportion_fibro.png", width = 6.2, height = 4, dpi = 300)
ggsave("celltype_proportion_fibro.pdf", width = 6.2, height = 4, dpi = 600)



nCells2 <- dplyr::count (Input.df, orig.ident)
nCells2 <- merge(nCells2, Input.df[,c(1,41)] %>% unique() %>% arrange(desc(chemo)), by = "orig.ident", all.x=TRUE)
pre_Counted2 <- Input.df %>% group_by(orig.ident) %>% dplyr::count(celltype)
Counted2 <- pre_Counted2 %>% ungroup %>% tidyr::complete(celltype, orig.ident,  fill = list(n = 0)) # with zero-values
Input2 <- merge (x=Counted2, y=nCells2, by = "orig.ident", all.x=TRUE)
Input2$Percent <- Input2$n.x / Input2$n.y * 100
Input2$chemo %>% table

my_comparisons <- list(c("Pre-treatment", "Post-treatment"))

Input2$chemo <- factor(Input2$chemo, levels = c("Pre-treatment", "Post-treatment"))       
ggplot(Input2, aes(chemo, Percent, fill = chemo)) + theme_bw() + NoLegend() + 
  geom_boxplot(width = 0.5) +
  geom_point(aes(color = chemo), size = 1.5, shape = 19) + geom_jitter(width = 0) + 
  labs(y = "% Proportions")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.text.x = element_text(angle = 50, vjust = 0.5),
        axis.text.y = element_text(size =15),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 15),   
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size=14)) + 
  facet_wrap(~ celltype, ncol = 3) +
  scale_fill_manual(values = chemo_color) + expand_limits(y=c(0, 105)) +
  scale_color_manual(values = chemo_color) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 5,label.y = 92)

ggsave("fibro_cluster_boxplot.png", width = 7, height = 6, dpi = 300)
ggsave("fibro_cluster_boxplot.pdf", width = 7, height = 6, dpi = 600)


saveRDS(fibro, "../fibroblast/OVARY5P_fibroblast_final.rds")


####omentum vs ovary 


fibro$metastatic_primary %>% table()
fibro$metastatic_primary <- fibro$orig.ident
new.meta.ids <- c('Metastatic','Primary','Primary','Primary','Primary','Primary','Primary','Primary','Metastatic','Primary')
fibro$metastatic_primary <- plyr::mapvalues(fibro$metastatic_primary, from = sample.list, to = new.meta.ids )

Idents(fibro) <- "orig.ident"
fibro_ov <- subset(fibro, idents = c("T110", "T92", "T17", "T18", "T34"))
fibro_om <- subset(fibro, idents = c("T04", "T70"))

fibro <- subset(fibro, idents = c("T110", "T92", "T17", "T18", "T34","T04", "T70"))
table(fibro$metastatic_primary, fibro$orig.ident)


Idents(fibro) <- "metastatic_primary"
clustering <- fibro@active.ident
cellsIn <- names(clustering[clustering == "Metastatic"])
cellsOut <- names(clustering[clustering == "Primary"])
fibro<- ScaleData(fibro)

fibro_markers <- FindMarkers(fibro, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)
fibro_markers$genes <- rownames(fibro_markers)

writexl::write_xlsx(fibro_markers, path = paste0("project_directory_name/ovary_project/Seurat/final/stromal/fibroblast/DEG_fibroblast_omentum_vs_ovary",".xlsx"))
fibro_markers$p_val_adj
EnhancedVolcano(fibro_markers, lab = rownames(fibro_markers), x = 'avg_log2FC', y = 'p_val_adj', title = "Omentum vs Ovary (Fibroblast)")

