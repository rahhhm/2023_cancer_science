library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

####Endothelial cells####
Idents(ovary) <- "celltype"

endo.path <- paste0(wd.path, "stromal/endothelial cells/")
setwd(endo.path)

endo <- subset(ovary, idents = "Endothelial cells")


table(endo$orig.ident, endo$chemo)


Idents(endo) <- "chemo"
clustering <- endo@active.ident
cellsIn <- names(clustering[clustering == "Post-treatment"])
cellsOut <- names(clustering[clustering == "Pre-treatment"])
endo<- ScaleData(endo)

endo_markers <- FindMarkers(endo, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0, )


endo_markers$significant <- ifelse(-log10(endo_markers$p_val_adj) > 5 & abs(endo_markers$avg_log2FC) >= 1,  
                                    ifelse(endo_markers$avg_log2FC >= 1 ,'Up' ,'Down'),' ')

writexl::write_xlsx(endo_markers, "OVARY5P_endo_treatment_DEG.xlsx" )

endo_markers$significant <- endo_markers$significant %>% as.factor()
endo_markers$significant <- factor(endo_markers$significant, levels = c("Up", "Down", " "))

ggplot(endo_markers, aes(avg_log2FC, -log10(p_val_adj), fill = significant, color = significant)) + 
  geom_point(size = 1.5, shape = 21) +
  theme_bw() + 
  geom_hline(yintercept=5, linetype='dashed', color='red', size=0.5) +
  geom_vline(xintercept=1, linetype = 'dashed', color='red', size = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color='red', size = 0.5) + ylim(c(0,20)) +
  scale_fill_manual(values=c("red","blue","gray67")) +
  scale_color_manual(values=c("red","blue","gray67")) +
  geom_text_repel(data = endo_markers[endo_markers[,"significant"] == "Up" | 
                                         endo_markers[,"significant"] == "Down" ,],
                  aes(label = genes), 
                  size = 5,
                  color = "black") + theme(legend.text = element_text(size = 15),
                                           legend.title = element_text(size = 15),
                                           axis.text = element_text(size = 15),
                                           axis.title = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA)))) 

ggsave("endo_markers_volcano.png", width = 8, height = 6, dpi = 300)
ggsave("endo_markers_volcano.pdf", width = 8, height = 6, dpi = 600)



endo <- NormalizeData(object=endo, normalization.method="LogNormalize", scale.factor=10000)

## Identification of highly variable features
endo <- FindVariableFeatures(endo, selection.method = "vst", nfeatures = num.features)
top10 <- head(VariableFeatures(endo), 10)

## Centering the data (z-scoring)
all.genes <- rownames(endo)
endo <- ScaleData(endo, features = all.genes)

## Run PCA
library(tibble)
endo <- RunPCA(endo, npcs=50, features = VariableFeatures(endo))

## Determine statistically significant principal components
#### Elbow plot
ElbowPlot(endo, ndims=50, reduction="pca") + theme_light()



##################################
#T107 excluded (cell = 1)
library(Seurat)
library(bbknnR)
#10, 0.8
Idents(endo) <- "orig.ident"
endo <- subset(endo, subset = orig.ident != "T107")
endo <- endo %>% 
  RunBBKNN(run_TSNE = FALSE, n_pcs = 22, batch_key = "orig.ident") %>% 
  FindClusters(resolution = 0.3, graph.name= "bbknn") 

DimPlot(endo, group.by = "chemo", label = TRUE, label.size = 5)

Idents(endo) <- "celltype"
endo_Markers = FindAllMarkers(object=endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
endo_Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
endo <- ScaleData(endo, features = top10$gene)


DoHeatmap(object=endo, features= top10$gene, draw.lines = FALSE, size = 4) + theme(text = element_text(size = 12)) + NoLegend()
ggsave(file="endo_marker_heatmap.png", width=6, height=8, dpi = 300) 


av_endo <-AverageExpression(endo, return.seurat = TRUE)

DoHeatmap(object=av_endo, features=top10$gene, draw.lines = FALSE, angle = 45, group.colors = celltype.colv,size = 4) + theme(text = element_text(size = 12)) + NoLegend()
ggsave(file="endo_average_heatmap.png", width=6.5, height=10, dpi = 300)
ggsave(file="endo_average_heatmap.pdf", width=6.5, height=10, dpi = 600)

celltype.id <- c("0"="EC0_ACKR1",  "1"="EC1_ANGPT2", "2"= "EC2_FABP4", 
                 "3"="EC3_ACTA2", "4"="EC4_PROX1")

identical(length(celltype.id), length(unique(endo$seurat_clusters)))
endo@meta.data$celltype <- plyr::mapvalues(endo@meta.data$seurat_clusters, from=names(celltype.id), to=celltype.id)
endo$celltype %>% table

pt.size <- 5

p1 <- DimPlot(endo, reduction = "umap", group.by = "celltype", pt.size= pt.size, cols = celltype.colv) + 
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
ggsave("dimplot_celltype_endo.png", width = 4, height = 4, dpi = 300)
ggsave("dimplot_celltype_endo.pdf", width = 4, height = 4, dpi = 600)



legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_endo_legend.png", width = 2, height = 2, dpi = 300)
ggsave("dimplot_celltype_endo_legend.pdf", width = 2, height = 2, dpi = 600)

p2 <- DimPlot(endo, reduction = "umap", group.by = "orig.ident", pt.size= pt.size, cols = orig.ident.colv) + 
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
ggsave("endo_dimplot_celltype_orig.png", width = 4, height = 4, dpi = 300)
ggsave("endo_dimplot_celltype_orig.pdf", width = 4, height = 4, dpi = 600)

legend <- get_legend(p2 + theme(legend.background = element_blank(),
                                legend.title = element_blank()) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend, xlim = c(1,1) )
ggsave("endo_dimplot_celltype_orig_legend.png", width = 1.5, height = 3.5, dpi = 300)
ggsave("endo_dimplot_celltype_orig_legend.pdf", width = 1.5, height = 3.5, dpi = 600)

endo$chemo <- factor(endo$chemo, levels =  c("Pre-treatment", "Post-treatment"))
p6 <- DimPlot(endo, reduction = "umap", group.by = "chemo", pt.size=pt.size, cols = c("#00B0F0", "#00B050")) + 
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
ggsave("endo_dimplot_total_treatment.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("../endothelial cells/endo_dimplot_total_treatment.pdf", width = 4, height = 4, dpi = 600)

legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("endo+dimplot_total_treatment_legend.png", width = 3, height = 1, dpi = 300)




Input.df <- endo@meta.data
Input.df <- Input.df %>% filter(celltype != "EC2_low-quality")
Input.df$celltype <- gdata::drop.levels(Input.df$celltype)

Input.df$orig.ident %>% table
Input.df$orig.ident <- factor(Input.df$orig.ident, levels = c("T70",
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

ggsave("seurat_clusters_proportion_endo.png", width = 6, height = 4, dpi = 300)
ggsave("seurat_clusters_proportion_endo.pdf", width = 6, height = 4, dpi = 600)

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
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(angle = 50, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 15),   
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size=14)) + 
  facet_wrap(~ celltype, ncol = 3) +
  scale_fill_manual(values = chemo_color) + expand_limits(y=c(0, 90)) +
  scale_color_manual(values = chemo_color) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 5,label.y = 80)

ggsave("endo_cluster_boxplot.png", width = 6, height = 5, dpi = 300)
ggsave("endo_cluster_boxplot.pdf", width = 6, height = 5, dpi = 600)


saveRDS(endo, "OVARY5P_endothelial_cell_final.rds")



####omentum vs ovary 


endo$metastatic_primary %>% table()
endo$metastatic_primary <- endo$orig.ident
new.meta.ids <- c('Metastatic','Primary','Primary','Primary','Primary','Primary','Primary','Primary','Metastatic','Primary')
endo$metastatic_primary <- plyr::mapvalues(endo$metastatic_primary, from = sample.list, to = new.meta.ids )

Idents(endo) <- "orig.ident"
endo_ov <- subset(endo, idents = c("T110", "T92", "T17", "T18", "T34"))
endo_om <- subset(endo, idents = c("T04", "T70"))

endo <- subset(endo, idents = c("T110", "T92", "T17", "T18", "T34","T04", "T70"))
table(endo$specimen_tissue_site, endo$orig.ident)


Idents(endo) <- "specimen_tissue_site"
clustering <- endo@active.ident
cellsIn <- names(clustering[clustering == "omentum"])
cellsOut <- names(clustering[clustering == "ovary"])
endo<- ScaleData(endo)

endo_markers <- FindMarkers(endo, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)

endo_markers$genes <- rownames(endo_markers)
writexl::write_xlsx(endo_markers, path = paste0("project_directory_name/ovary_project/Seurat/final/stromal/endothelial cells/DEG_endothelial_omentum_vs_ovary",".xlsx"))
endo_markers$p_val_adj
EnhancedVolcano(endo_markers, lab = rownames(endo_markers), x = 'avg_log2FC', y = 'p_val_adj', title = "Omentum vs Ovary (Endothelial cells)")


