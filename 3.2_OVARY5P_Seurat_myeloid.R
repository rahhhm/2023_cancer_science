####Myeloid cells####
library(RColorBrewer)
library(ggplot2)
library(Seurat)


myeloid.wd <- paste0(wd.path, "stromal/myeloid_cells")
dir.create(myeloid.wd)
setwd(myeloid.wd)

Idents(ovary_total) <- "celltype"
myeloid <- subset(ovary_total, idents = "Myeloid cells")

table(myeloid$orig.ident, myeloid$chemo)

Idents(myeloid) <- "chemo"
clustering <- myeloid@active.ident
cellsIn <- names(clustering[clustering == "Post-treatment"])
cellsOut <- names(clustering[clustering == "Pre-treatment"])
myeloid<- ScaleData(myeloid)

myeloid_markers <- FindMarkers(myeloid, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0, )

myeloid_markers$significant <- ifelse(-log10(myeloid_markers$p_val_adj) > 5 & abs(myeloid_markers$avg_log2FC) >= 1, 
                                      ifelse(myeloid_markers$avg_log2FC >= 1 , 'Up' ,'Down'),' ')

writexl::write_xlsx(myeloid_markers, "OVARY5P_myeloid_treatment_DEG.xlsx" )


myeloid_markers$significant <- myeloid_markers$significant %>% as.factor()
myeloid_markers$significant <- factor(myeloid_markers$significant, levels = c("Up", "Down", " "))

ggplot(myeloid_markers, aes(avg_log2FC, -log10(p_val_adj), fill = significant, color = significant)) + 
  geom_point(size = 1.5, shape = 21) +
  theme_bw() + 
  geom_hline(yintercept=5, linetype='dashed', color='red', size=0.5) +
  geom_vline(xintercept=1, linetype = 'dashed', color='red', size = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color='red', size = 0.5) + ylim(c(0,300)) +
  scale_fill_manual(values=c("red","blue","gray67")) +
  scale_color_manual(values=c("red","blue","gray67")) +
  geom_text_repel(data = myeloid_markers[myeloid_markers[,"significant"] == "Up" | 
                                           myeloid_markers[,"significant"] == "Down" ,],
                  aes(label = genes), 
                  size = 5,
                  color = "black") + theme(legend.text = element_text(size = 15),
                                           legend.title = element_text(size = 15),
                                           axis.text = element_text(size = 15),
                                           axis.title = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA)))) 
ggsave("myeloid_markers_volcano.png", width = 8, height = 6, dpi = 300)
ggsave("myeloid_markers_volcano.pdf", width = 8, height = 6, dpi = 600)



myeloid <- NormalizeData(object=myeloid, normalization.method="LogNormalize", scale.factor=10000)

## Identification of highly variable features
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = num.features)
top10 <- head(VariableFeatures(myeloid), 10)

## Centering the data (z-scoring)
all.genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid)

## Run PCA
library(tibble)
myeloid <- RunPCA(myeloid, npcs=50, features = VariableFeatures(myeloid))

## Determine statistically significant principal components
#### Elbow plot
ElbowPlot(myeloid, ndims=50, reduction="pca") + theme_light()

##################################
library(Seurat)
library(bbknnR)

myeloid <- myeloid %>% 
  RunBBKNN(run_TSNE = FALSE, n_pcs = 20, batch_key = "orig.ident") %>% 
  FindClusters(resolution = 0.3, graph.name= "bbknn") 

Idents(myeloid) <- "celltype"
myeloid_Markers = FindAllMarkers(object=myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
myeloid_Markers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) -> top10
myeloid <- ScaleData(myeloid, features = top10$gene)


DoHeatmap(object=myeloid, features= top10$gene, draw.lines = FALSE, size = 4) + theme(text = element_text(size = 12)) + NoLegend()
ggsave(file="endo_marker_heatmap.png", width=6, height=8, dpi = 300) 




av_myeloid <-AverageExpression(myeloid, return.seurat = TRUE)

DoHeatmap(object=av_myeloid, features=top10$gene, draw.lines = FALSE, angle = 60, group.colors = celltype.colv,size = 4) + theme(text = element_text(size = 12)) + NoLegend()
ggsave(file="myeloid_average_heatmap.png", width=8, height=10, dpi = 300)
ggsave(file="myeloid_average_heatmap.pdf", width=8, height=10, dpi = 600)


myeloid$nFeature_RNA
VlnPlot(myeloid, features = c("nCount_RNA", "nFeature_RNA"))

celltype.id <- c("0"="MC0_APOE",  "3"="MC1_STMN1",
                 "2"= "MC2_FCN1", "5"="MC3_LAMP3", "6"="MC4_IRF4", "4"="MC5_Doublet","1"="MC6_Low-quality")

identical(length(celltype.id), length(unique(myeloid$seurat_clusters)))
myeloid@meta.data$celltype <- plyr::mapvalues(myeloid@meta.data$seurat_clusters, from=names(celltype.id), to=celltype.id)
myeloid$celltype %>% table

myeloid$celltype <- factor(myeloid$celltype, levels = c("MC0_APOE", "MC1_STMN1",
                                                        "MC2_FCN1", "MC3_LAMP3",
                                                        "MC4_IRF4","MC5_Doublet","MC6_Low-quality"))

Idents(myeloid) <- "celltype"
pt.size <- 2
p1 <- DimPlot(myeloid, reduction = "umap", group.by = "celltype", pt.size= pt.size, cols = celltype.colv) + 
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

ggsave("dimplot_celltype_myeloid.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("dimplot_celltype_myeloid.pdf", width = 5.5, height = 5.5, dpi = 600)

legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_myeloid_legend.png", width = 2, height = 3, dpi = 300)
ggsave("dimplot_celltype_myeloid_legend.pdf", width = 2, height = 3, dpi = 600)

p2 <- DimPlot(myeloid, reduction = "umap", group.by = "orig.ident", pt.size= pt.size, cols = orig.ident.colv) + 
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
ggsave("myeloid_dimplot_celltype_orig.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("myeloid_dimplot_celltype_orig.pdf", width = 5.5, height = 5.5, dpi = 600)

legend <- get_legend(p2 + theme(legend.background = element_blank(),
                                legend.title = element_blank()) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend, xlim = c(1,1) )
ggsave("myeloid_dimplot_celltype_orig_legend.png", width = 1.5, height = 3.5, dpi = 300)
ggsave("myeloid_dimplot_celltype_orig_legend.pdf", width = 1.5, height = 3.5, dpi = 600)
myeloid$chemo <- factor(myeloid$chemo, levels = c("Pre-treatment", "Post-treatment"))
p6 <- DimPlot(myeloid, reduction = "umap", group.by = "chemo", pt.size=pt.size, cols = c("#00B0F0", "#00B050")) + 
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
ggsave("myeloid_dimplot_total_treatment.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("../myeloid_cells/myeloid_dimplot_total_treatment.pdf", width = 5.5, height = 5.5, dpi = 600)


legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("myeloid+dimplot_total_treatment_legend.png", width = 3, height = 1, dpi = 300)




Input.df <- myeloid@meta.data
Input.df <- Input.df %>% filter(celltype != "MC5_Doublet" & celltype != "MC6_Low-quality")
Input.df$celltype <- gdata::drop.levels(Input.df$celltype)

Input.df$orig.ident %>% table
Input.df$orig.ident <- factor(Input.df$orig.ident, levels = c("T107", "T70",
                                                              "T04", "T110",
                                                              "T92", "T66",
                                                              "T17", "T18",
                                                              "T39", "T34"))
nCells2 <- dplyr::count (Input.df, orig.ident)
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

ggsave("seurat_clusters_proportion_myeloid.png", width = 6, height = 4, dpi = 300)
ggsave("seurat_clusters_proportion_myeloid.pdf", width = 6, height = 4, dpi = 600)

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
        axis.text.y = element_text(size =12),
        axis.text.x = element_text(angle = 50, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 12),   
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size=14)) + 
  facet_wrap(~ celltype, ncol = 3) +
  scale_fill_manual(values = chemo_color) + expand_limits(y=c(0, 95)) +
  scale_color_manual(values = chemo_color) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 5,label.y = 85)

ggsave("myeloid_cluster_boxplot.png", width = 5, height = 4.5, dpi = 300)
ggsave("myeloid_cluster_boxplot.pdf", width = 5, height = 4.5, dpi = 600)

saveRDS(myeloid, "OVARY5P_final_myeloid.rds")

Idents(myeloid)

#mac, OV V OMENTUM

mac <- subset(myeloid, idents = "MC0_APOE")
dim(mac)

mac$orig.ident %>% table

mac$metastatic_primary %>% table()
mac$metastatic_primary <- mac$orig.ident
new.meta.ids <- c('Metastatic','Primary','Primary','Primary','Primary','Primary','Primary','Primary','Metastatic','Primary')
mac$metastatic_primary <- plyr::mapvalues(mac$metastatic_primary, from = sample.list, to = new.meta.ids )

Idents(mac) <- "orig.ident"
mac_ov <- subset(mac, idents = c("T110", "T92", "T17", "T18", "T34"))
mac_om <- subset(mac, idents = c("T04", "T70"))

mac <- subset(mac, idents = c("T110", "T92", "T17", "T18", "T34","T04", "T70"))
table(mac$specimen_tissue_site, mac$orig.ident)


Idents(mac) <- "metastatic_primary"
clustering <- mac@active.ident
cellsIn <- names(clustering[clustering == "Metastatic"])
cellsOut <- names(clustering[clustering == "Primary"])
mac<- ScaleData(mac)

mac_markers <- FindMarkers(mac, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)
mac_markers$genes <- rownames(mac_markers)

writexl::write_xlsx(mac_markers, path = paste0("project_directory_name/ovary_project/Seurat/final/stromal/myeloid_cells/DEG_macrophage_omentum_vs_ovary",".xlsx"))

mac_markers$p_val_adj
EnhancedVolcano(mac_markers, lab = rownames(mac_markers), x = 'avg_log2FC', y = 'p_val_adj', title = "Omentum vs Ovary")

Idents(myeloid) <- "orig.ident"
table(myeloid$specimen_tissue_site, myeloid$orig.ident)

myeloid_om <- subset(myeloid, idents = c("T110", "T92", "T17", "T18", "T34","T04", "T70"))
Idents(myeloid_om) <- "specimen_tissue_site"
clustering <- myeloid_om@active.ident
cellsIn <- names(clustering[clustering == "omentum"])
cellsOut <- names(clustering[clustering == "ovary"])
myeloid_om<- ScaleData(myeloid_om)

myeloid_om_markers <- FindMarkers(myeloid_om, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)
myeloid_om_markers$genes <- rownames(myeloid_om_markers)
writexl::write_xlsx(myeloid_om_markers, path = paste0("project_directory_name/ovary_project/Seurat/final/stromal/myeloid_cells/DEG_myeloid_omentum_vs_ovary",".xlsx"))
EnhancedVolcano(myeloid_om_markers, lab = rownames(myeloid_om_markers), x = 'avg_log2FC', y = 'p_val_adj', title = "Omentum vs Ovary")

ovary



