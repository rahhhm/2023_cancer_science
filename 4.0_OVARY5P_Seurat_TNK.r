############T/NK subclustering###############

####LOADING total cell######
library(Seurat)
library(ggplot2)
library(bbknnR)
library(dplyr)
                  
label = "OVARY5P"
wd.path = "D:/project_directory_name/ovary_project/Seurat/final/"
tnk.path = paste0(wd.path, "tnk/")
dir.create(tnk.path)
setwd(tnk.path)
getwd()
####T/NK cells####
Idents(ovary) <- "celltype"
tnk <- subset(ovary, idents = "T/NK cells")

####Treatment DEG####
Idents(tnk) <- "chemo"
clustering <- tnk@active.ident
cellsIn <- names(clustering[clustering == "Post-treatment"])
cellsOut <- names(clustering[clustering == "Pre-treatment"])
tnk<- ScaleData(tnk)

tnk_markers <- FindMarkers(tnk, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)

tnk_markers$significant <- ifelse(-log10(tnk_markers$p_val_adj) > 100 & abs(tnk_markers$avg_log2FC) >= 1, 
                                  ifelse(tnk_markers$avg_log2FC >= 1 , 'Up' ,'Down'),' ')


tnk_markers$significant <- tnk_markers$significant %>% as.factor()
tnk_markers$significant <- factor(tnk_markers$significant, levels = c("Up", "Down", " "))
writexl::write_xlsx(tnk_markers, "OVARY5P_tnk_treatment_DEG.xlsx" )

ggplot(tnk_markers, aes(avg_log2FC, -log10(p_val_adj), fill = significant, color = significant)) + 
  geom_point(size = 1.5, shape = 21) +
  theme_bw() + 
  geom_hline(yintercept=100, linetype='dashed', color='red', size=0.5) +
  geom_vline(xintercept=1, linetype = 'dashed', color='red', size = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color='red', size = 0.5) + ylim(c(0,300)) +
  scale_fill_manual(values=c("red","blue","gray67")) +
  scale_color_manual(values=c("red","blue","gray67")) +
  geom_text_repel(data = tnk_markers[tnk_markers[,"significant"] == "Up" | 
                                       tnk_markers[,"significant"] == "Down" ,],
                  aes(label = genes), 
                  size = 5,
                  color = "black") + theme(legend.text = element_text(size = 15),
                                           legend.title = element_text(size = 15),
                                           axis.text = element_text(size = 15),
                                           axis.title = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA)))) 
ggsave("tnk_markers_volcano.png", width = 8, height = 6, dpi = 300)
ggsave("tnk_markers_volcano.pdf", width = 8, height = 6, dpi = 600)


#################################################################################################################

Idents(ovary) <- "celltype"
tnk <- subset(ovary, idents = "T/NK cells")
####T/NK cells subset####
# for variable genes
num.features = 3000

## Log Normalized data (TPM-like values)
as.numeric(tnk@assays$RNA@counts) %>%  summary

tnk@assays$RNA@counts %>%  head(20)

tnk <- NormalizeData(object=tnk, normalization.method="LogNormalize", scale.factor=10000)

## Identification of highly variable features
tnk <- FindVariableFeatures(tnk, selection.method = "vst", nfeatures = num.features)
top10 <- head(VariableFeatures(tnk), 10)

## Centering the data (z-scoring)
all.genes <- rownames(tnk)
tnk <- ScaleData(tnk, features = all.genes, vars.to.regress = "percent.mt")

## Run PCA
library(tibble)
tnk <- RunPCA(tnk, npcs=50, features = VariableFeatures(tnk))

  ## Determine statistically significant principal components
#### Elbow plot
ElbowPlot(tnk, ndims=50, reduction="pca") + theme_light()

##################################
library(Seurat)
library(bbknnR)

tnk3 <- tnk %>% 
  RunBBKNN(run_TSNE = FALSE, n_pcs = 13, batch_key = "orig.ident") %>% 
  FindClusters(resolution = 1.5, graph.name= "bbknn") 

################DEG###################
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gghighlight)
library(tidyr)
library(readxl)

## Identify DEG
getwd()
Idents(tnk) <- "seurat_clusters"
all_Markers = FindAllMarkers(object=tnk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

writexl::write_xlsx(all_Markers, path = paste0("cluster_markers_Ttest_", "_pc", 13, "_res", 1.5, ".xlsx"))
all_Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

writexl::write_xlsx(top10, path = paste0("cluster_markers_Ttest_", "_pc", 13, "_res", 1.5,"top10", ".xlsx"))

tnk <- ScaleData(tnk, features = top10$gene)
DoHeatmap(object=tnk, features=top10$gene) + NoLegend() + theme(text = element_text(size = 10))
OutPutFile_Marker_heatmap = paste0("Marker_heatmap_Ttest_", ".png")
ggsave(file=OutPutFile_Marker_heatmap, width=33, height=15) 
getwd()

###average_heatmap###
Idents(tnk) <- "celltype"
tnk_Markers = FindAllMarkers(object=tnk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
tnk_Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
tnk <- ScaleData(tnk, features = avg_feature)

av_tnk <-AverageExpression(tnk, return.seurat = TRUE)

DoHeatmap(object= av_tnk, features = top10$gene, draw.lines = FALSE, size = 5) + theme(text = element_text(size = 13)) + NoLegend()

ggsave(file="tnk_average_heatmap.png", width= 20, height= 15, dpi = 300)

##### 5. Define cell types ###############################################################################
## Write your cell type name for each cluster
library(RColorBrewer)

celltype.id <- c("0"="CD8 Tex", "1"="NK", "2"= "CD4 Tm", "3"="CD8 Tem",
                 "4"="CD8 T ISG", "5"="CD4 Treg", "6"="CD4 Tn", "7"="CD4 Tfh",
                 "8"="Proliferating", "9"="CD4 Treg", "10"="CD8 Tem", "11"="CD4 Tm", 
                 "12" = "CD8 Tn/m", "13" = "CD4 T ISG", "14" = "NK")

identical(length(celltype.id), length(unique(tnk$seurat_clusters)))
tnk@meta.data$celltype <- plyr::mapvalues(tnk@meta.data$seurat_clusters, from=names(celltype.id), to=celltype.id)
tnk$celltype %>% table

tnk$celltype <- factor(tnk$celltype, levels = c("CD4 Tn", "CD4 Tm", "CD4 T ISG", "CD4 Tfh",
                                                  "CD4 Treg", "CD8 Tn/m", "CD8 Tem",
                                                  "CD8 T ISG", "CD8 Tex", "NK", "Proliferating"))

tnk$orig.ident <- factor(tnk$orig.ident, levels = rev(c("T107", "T70",
                                                          "T04", "T110",
                                                          "T92", "T66",
                                                          "T17", "T18",
                                                          "T39", "T34")))

celltype.colv <- c("CD4 Tn" = "#7ece00",
                    "CD4 Tm" = "#00b10f",
                    "CD4 T ISG" = "#e466b3",
                    "CD4 Tfh" = "#bc6c25",
                    "CD4 Treg" = "#6d2727",
                    "CD8 Tn/m" = "#ff882c",
                    "CD8 Tem" = "#0BC7D2",
                    "CD8 T ISG" = "#AE64C6",
                    "CD8 Tex" = "#e92244",
                    "NK" = "#27187e", 
                    "Proliferating" = "#2a9d8f")

orig.ident.colv <- c("T107" = "#5C5EAC",
                     "T70" = "#5377AE", 
                     "T04" = "#114292",
                     "T110" = "#A779C5",
                     "T92" = "#FA6D3E",  
                     "T66" = "#E34F79", 
                     "T17" = "#FE9920", 
                     "T18" = "#F4CB2A", 
                     "T39" = "#0BC7D2", 
                     "T34" = "#ABD84C")


celltype_table <- table(tnk$orig.ident, tnk$celltype) %>%  as.matrix.data.frame()
row.names(celltype_table) <- names(tnk$orig.ident%>% table)
colnames(celltype_table) <- names(tnk$celltype%>% table)
celltype_table <- celltype_table[,-11]
write.table(celltype_table, "celltype_table.csv", sep = ",", row.names = T, col.names = NA)

##########################################################################################

pt.size <- 1

p1 <- DimPlot(tnk, reduction = "umap", group.by = "seurat_clusters", pt.size= pt.size) + 
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
ggsave("dimplot_seurat_clusters_tnk.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("dimplot_seurat_clusters_tnk.pdf", width = 5.5, height = 5.5, dpi = 600)


legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_seurat_clusters_tnk_legend.png", width = 2, height = 3, dpi = 300)
ggsave("dimplot_seurat_clusters_tnk_legend.pdf", width = 2, height = 3, dpi = 600)

p2 <- DimPlot(tnk, reduction = "umap", group.by = "orig.ident", pt.size= pt.size, cols = orig.ident.colv) + 
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
ggsave("tnk_dimplot_celltype_orig.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("tnk_dimplot_celltype_orig.pdf", width = 5.5, height = 5.5, dpi = 600)

legend <- get_legend(p2 + theme(legend.background = element_blank(),
                                legend.title = element_blank()) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend, xlim = c(1,1) )
ggsave("tnk_dimplot_celltype_orig_legend.png", width = 1.5, height = 3.5, dpi = 300)

p6 <- DimPlot(tnk, reduction = "umap", group.by = "chemo", pt.size=pt.size, cols = c("#00B0F0", "#00B050")) + 
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
ggsave("tnk_dimplot_total_treatment.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("tnk_dimplot_total_treatment.pdf", width = 5.5, height = 5.5, dpi = 600)



legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("tnk+dimplot_total_treatment_legend.png", width = 3, height = 1, dpi = 300)


p1 <- DimPlot(tnk, reduction = "umap", group.by = "celltype", pt.size= pt.size, cols = celltype.colv) + 
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
ggsave("dimplot_celltype_tnk.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("dimplot_celltype_tnk.pdf", width = 5.5, height = 5.5, dpi = 600)

legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_tnk_legend.png", width = 2, height = 3.5, dpi = 300)
ggsave("dimplot_celltype_tnk_legend.pdf", width = 2, height = 3.5, dpi = 600)

library(pheatmap)

Idents(tnk) <- "seurat_clusters"

tnk$celltype <- factor(tnk$celltype, levels = c("CD4 Tn", "CD4 Tm", "CD4 T ISG","CD4 Tfh", "CD4 Treg","CD8 Tn/m", "CD8 Tem",  "CD8 T ISG", "CD8 Tex", "NK", "Proliferating"))


tnk$celltype %>% table
Idents(tnk) <- "celltype"

target.list<- list(CD4_8_features <- c("CD4", "CD40LG","CD8A", "CD8B"),
naive_features <- c("CCR7", "SELL"),
memory_features  <- c("IL7R", "GPR183", "LMNA") ,
ISG_features <- c("IFIT3", "ISG15"),
follicular_helper_features <- c("CD200", "CXCL13"),
regulatory_features <- c("TNFRSF9", "CTLA4", "FOXP3", "IL2RA", "TIGIT") ,
effector_features <- c( "KLRG1", "GZMK", "GZMA", "GZMH"),
exhausted_features <- c("HAVCR2", "GZMB" , "LAG3", "PDCD1"),
NK_features <-c("FGFBP2", "GNLY", "NCAM1", "XCL1"),
proliferating_features <-c("MKI67", "TOP2A"))

n_height<- c(2,4,4,4,5,2,2,3,2,4)

#+ scale_size(breaks = c(0, 25, 50)) 
d1 <- DotPlot(tnk, features = proliferating_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend() + theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold")) + labs(x = paste0("Proliferating"))

d2 <- DotPlot(tnk, features = NK_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend() + theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold")) + labs(x = paste0("NK"))

d3 <- DotPlot(tnk, features = exhausted_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold"))+ labs(x = paste0("Exhausted"))

d4 <- DotPlot(tnk, features = effector_mermory_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend()+ theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold"))+ labs(x = "Effector")


d5 <- DotPlot(tnk, features = regulatory_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend()+ theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold"))+ labs(x = paste0("Regulatory"))

d6 <-  DotPlot(tnk, features = follicular_helper_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend()+ theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold"))+ labs(x = paste0("Follicular", '\n',"helper"))

d7 <-DotPlot(tnk, features = ISG_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend()+ theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(12), title = element_text(size = 11, face = "bold"))+ labs(x = "ISG")

d8 <-  DotPlot(tnk, features = memory_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend()+ theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold"))+ labs(x = "Memory")

d9 <-  DotPlot(tnk, features = naive_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend()+ theme(axis.line.x = element_blank(),axis.title.x = element_blank(), axis.title.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 12), title = element_text(size = 11, face = "bold")) + labs(x = "Naive")

d0 <- DotPlot(tnk, features = CD4_8_features, split.by = "seurat_clusters", cols = "RdYlBu") + coord_flip() + NoLegend() +  
  labs(x = "CD4_8") + theme(axis.text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title.x.bottom = element_blank(), title = element_text(size = 11, face = "bold")) 

d1 + d2 + d3 + d4 + d5 +  d6 + d7 + d8 + d9 + d0 + patchwork::plot_layout(ncol = 1, heights = n_height) 

ggsave("dotplot_celltype_tnk.png", width = 8, height = 13.5, dpi = 300)
DimPlot(tnk, group.by = "seurat_clusters" )

tnk$cluster_celltype <- paste0(tnk$celltype, "_", tnk$seurat_clusters)
Idents(tnk) <- "cluster_celltype"
tnk$cluster_celltype <- factor(tnk$cluster_celltype, 
                                       levels = rev(c("CD4 Tn_6",
                                                      "CD4 Tm_2", "CD4 Tm_11",
                                                      "CD4 T ISG_13", "CD4 Tfh_7","CD4 Treg_5", "CD4 Treg_9",
                                                      "CD8 Tn/m_12", "CD8 Tem_3", "CD8 Tem_10", "CD8 T ISG_4",
                                                      "CD8 Tex_0", "NK_1", "NK_14", "Proliferating_8")))

DotPlot(tnk, features = rev(paste(unlist(target.list))), 
        cols= "RdYlBu", 
        dot.min = 0, dot.scale = 5.5) + RotatedAxis() +
  scale_x_discrete(position='top') +
  scale_y_discrete(position='right') +
  theme_bw() + theme(legend.position='right',
                     legend.direction = "vertical",
                     legend.key.width = unit(0.22,"cm"),
                     legend.key.size =  unit(0.2,"cm"),
                     legend.background = element_blank(),
                     legend.title = element_text(colour="black", size=11),
                     legend.text = element_text(colour="black", size=10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.title = element_blank(),
                     axis.text.x = element_text(colour="black", size=11, angle=45, hjust=0, vjust=1),
                     axis.text.y = element_text(colour="black", size=11))

ggsave("dotplot_celltype_tnk.png", width = 10, height = 4.5, dpi = 300)
ggsave("D:/project_directory_name/ovary_project/Seurat/final/tnk/dotplot_celltype_tnk.pdf", width = 10, height = 4.5, dpi = 600)

p5 <-  FeaturePlot(tnk, features = "ZNF683", pt.size = 2) +
  theme(aspect.ratio = 1,
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

AugmentPlot(p5)
ggsave("featureplot_ZNF683.png", width = 3, height = 3, dpi = 300)
ggsave("../featureplot_ZNF683.pdf", width = 3, height = 3, dpi = 600)

legend <- get_legend(p5)

ggdraw(legend, xlim = c(1,1))
ggsave("../featureplot_GNLY_legend.png", width = 2, height = 2, dpi = 300)
ggsave("../featureplot_ZNF683_legend.pdf", width = 2, height = 2, dpi = 600)


p0 <- FeaturePlot(tnk, features = c("IFIT3", "GZMK", "GZMB","LAG3", "PRF1", "GNLY"), ncol = 3, pt.size = pt.size) &
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
AugmentPlot(p0)


ggsave("fig3e.png", width = 6.75, height = 4, dpi = 300)
ggsave("fig3e.pdf", width = 6.75, height = 4, dpi = 300)



legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_tnk_legend.png", width = 2, height = 3.5, dpi = 300)




p5 <-  FeaturePlot(tnk, features = "TOX", pt.size = 2) +
  theme(aspect.ratio = 1,
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

AugmentPlot(p5)
ggsave("featureplot TOX .pdf", width = 3, height = 3, dpi = 600)

legend <- get_legend(p5)

ggdraw(legend, xlim = c(1,1))
ggsave("featureplot TOX legend.pdf", width = 2, height = 2, dpi = 600)


c("PDCD1", "HAVCR2", "TIGIT", "TOX")





P1 <- FeaturePlot(tnk, features = c("CD3E","MKI67", "TRDC", "NCAM1"), ncol = 4) & labs(title = "")  & 
  theme(aspect.ratio = 1,
        legend.position='none', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10)) 

AugmentPlot(P1)

ggsave("featureplot_marker1.png", width = 8, height = 2, dpi = 300)
ggsave("featureplot_marker1.pdf", width = 8, height = 2, dpi = 600)


FeaturePlot(tnk, features = c("SELL", "CCR7", "IL7R", 'GPR183'), ncol = 4) &
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10))

ggsave("featureplot_marker2.png", width = 8, height = 2, dpi = 300)


FeaturePlot(tnk, features = c('CXCL13', 'CD200', "FOXP3", "TNFRSF9"), ncol = 4) &
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10))

ggsave("featureplot_marker3.png", width = 8, height = 2, dpi = 300)

FeaturePlot(tnk, features = c("KLRG1", "GZMH", "GZMA", "NKG7"), ncol = 4) &
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10))

ggsave("featureplot_marker4.png", width = 8, height = 2, dpi = 300)

FeaturePlot(tnk, features = c("PDCD1", "HAVCR2", "TIGIT", "TOX"), ncol = 4) &
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10))

ggsave("featureplot_marker5.png", width = 8, height = 2, dpi = 300)







FeaturePlot(tnk, features=c('IFIT3', 'ISG15', "CXCL13", "CD200", "FOXP3", "TNFRSF9"), pt.size=0.01, ncol = 4) &
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10))
ggsave("featureplot_marker_CD4.png", width = 9, height = 4, dpi = 300)

FeaturePlot(tnk, features = c("PDCD1","TIGIT","CTLA4"), ncol = 3) &
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10))
ggsave("D:/project directory name/Seurat/final/BBKNN/paper_figure/featureplot_marker_CD4.png", width = 6.75, height = 2, dpi = 300)



P0 <- FeaturePlot(tnk, features = c("ITGAE", "ZNF683"), ncol = 2) &
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10)) 

brewer.pal(n = 2, name = "Dark2")
Idents(tnk) <- "celltype"
VlnPlot(tnk, features = c("ITGAE", "ZNF683"), stack = T, flip = T, idents = c("CD4 Tn", "CD4 T ISG", "CD8 Tn/m", "CD4 Tm", "CD4 Treg", "CD8 Tem", "CD8 T ISG", "CD4 Tfh","Proliferating","CD8 Tex"), cols = c("#1B9E77", "#D95F02")) + NoLegend() +
  theme(axis.title.x.bottom = element_blank(),
        axis.text.x = element_text(size = 8))

ggsave("D:/project directory name/Seurat/final/BBKNN/paper_figure/trm_marker.png", width = 9, height = 3, dpi = 300)
ggsave("../trm_marker.pdf", width = 6, height = 3, dpi = 600)


FeaturePlot(tnk, features = c("CXCR6")) +
  theme(aspect.ratio = 1,
        legend.position='right', 
        legend.key.size = unit(0.4, 'cm'),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        title = element_text(size = 10)) |
  VlnPlot(tnk, features = c("CXCR6"), flip = T, pt.size = 0) + NoLegend() +
  theme(axis.title.x.bottom = element_blank(),
        axis.text.x = element_text(size = 8))

ggsave("D:/project directory name/Seurat/final/BBKNN/paper_figure/cxcr6_marker.png", width = 6, height = 3, dpi = 300)

# Input.df <- tnk@meta.data
# Input.df %>% colnames
# nCells2 <- dplyr::count (Input.df, orig.ident)
# nCells2 <- merge(nCells2, Input.df[,c(1,20)] %>% unique() %>% arrange(desc(chemo)), by = "orig.ident", all.x=TRUE)
# pre_Counted2 <- Input.df %>% group_by(orig.ident) %>% dplyr::count(celltype)
# Counted2 <- pre_Counted2 %>% ungroup %>% tidyr::complete(celltype, orig.ident,  fill = list(n = 0)) # with zero-values
# Input2 <- merge (x=Counted2, y=nCells2, by = "orig.ident", all.x=TRUE)
# Input2$Percent <- Input2$n.x / Input2$n.y * 100
# Input2$chemo %>% table
# 
# Input3 <- Input2 %>% filter(orig.ident != "T39")
# Input4 <- Input3 %>% filter(orig.ident != "T107")
# 
# my_comparisons <- list(c("Pre-treatment", "Post-treatment"))
# 
# Input4$chemo <- factor(Input4$chemo, levels = c("Pre-treatment", "Post-treatment"))       
# ggplot(Input4, aes(chemo, Percent, fill = chemo)) + theme_bw() + NoLegend() + 
#   geom_boxplot(width = 0.5) +
#   geom_point(aes(color = chemo), size = 1.5, shape = 19) + geom_jitter(width = 0) + 
#   labs(y = "% Proportions")+
#   theme(panel.background = element_rect(inherit.blank = F),
#         axis.text.x = element_text(angle = 50, vjust = 0.5),
#         axis.ticks.x = element_blank(),
#         axis.title.x.bottom = element_blank(),
#         axis.title.y.left = element_text(size = 12),   
#         axis.text.x.bottom = element_blank(),
#         strip.text = element_text(size=14)) + 
#   facet_wrap(~ celltype, ncol = 4) +
#   scale_fill_manual(values = chemo_color) + expand_limits(y=c(0, 75)) +
#   scale_color_manual(values = chemo_color) + 
#   stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4,label.y = 70)
# 
# markers <- list()
# markers$cytotoxic <- c("PRF1", "GZMH", "GNLY", "IFNG")
# markers$dysfunction <- c("LAG3", "PDCD1", "TIGIT", "GZMB")
# markers$naive <- c("CCR7", "IL7R", "GPR183", "LMNA" )
# markers$regulatory <- c("IL2RA","TNFRSF9", "FOXP3", "CTLA4")
# markers$ifn <- c("IFI27", "MX1", "OAS2", "ISG15", "IFIT3", "STAT1")
# names(markers) <- c("cytotoxic", "exhausted", "naive", "regulatory", "ifn")

Idents(tnk) <- "orig.ident"
table(tnk$specimen_tissue_site, tnk$orig.ident)

tnk_om <- subset(tnk, idents = c("T110", "T92", "T17", "T18", "T34","T04", "T70"))
Idents(tnk_om) <- "specimen_tissue_site"
clustering <- tnk_om@active.ident
cellsIn <- names(clustering[clustering == "omentum"])
cellsOut <- names(clustering[clustering == "ovary"])
tnk_om<- ScaleData(tnk_om)
tnk_om$orig.ident %>% table
tnk_om_markers <- FindMarkers(tnk_om, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)
tnk_om_markers$genes <- rownames(tnk_om_markers)
writexl::write_xlsx(tnk_om_markers, path = paste0("project_directory_name/ovary_project/Seurat/final/stromal/tnk_cells/DEG_tnk_omentum_vs_ovary",".xlsx"))
EnhancedVolcano(tnk_om_markers, lab = rownames(tnk_om_markers), x = 'avg_log2FC', y = 'p_val_adj', title = "Omentum vs Ovary")
ggsave(paste0("T110 vs T92" ,".png"), width = 15, height = 12)


####ISG15 scoring####
library(cowplot)

tnk <- AddModuleScore(tnk, features = list(c("IFIT3", "ISG15", "MX1")))

p5 <- FeaturePlot(tnk, "Cluster1", min.cutoff = 0, pt.size = 2) + labs(title = "ISG score") +
  theme(aspect.ratio = 1,
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

AugmentPlot(p5)
ggsave("d:/project_directory_name/ovary_project/Seurat/final/tnk/featureplot_ISG_Score.png", width = 3, height = 3, dpi = 300)
ggsave("d:/project_directory_name/ovary_project/Seurat/final/tnk/featureplot_ISG_Score.pdf", width = 3, height = 3, dpi = 600)

legend <- get_legend(p5)

ggdraw(legend, xlim = c(1,1))
ggsave("d:/project_directory_name/ovary_project/Seurat/final/tnk/featureplot_ISG_Score_legend.pdf", width = 2, height = 2, dpi = 600)


