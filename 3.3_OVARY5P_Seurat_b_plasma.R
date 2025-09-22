
####B/Plasma cells####
b_plasma.wd <- paste0(wd.path, "stromal/b_plasma_cells")
dir.create(b_plasma.wd)
setwd(b_plasma.wd)

Idents(ovary_total) <- "celltype"
b_plasma <- subset(ovary_total, idents = "B/Plasma cells")

table(b_plasma$orig.ident, b_plasma$chemo)
Idents(b_plasma) <- "orig.ident"

b_plasma_e107 <- subset(b_plasma, subset = orig.ident != "T107")
table(b_plasma_e107$orig.ident, b_plasma_e107$chemo)

Idents(b_plasma_e107) <- "chemo"
clustering <- b_plasma_e107@active.ident
cellsIn <- names(clustering[clustering == "Post-treatment"])
cellsOut <- names(clustering[clustering == "Pre-treatment"])
b_plasma_e107<- ScaleData(b_plasma_e107)

b_plasma_e107_markers <- FindMarkers(b_plasma_e107, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)

b_plasma_e107_markers$significant <- ifelse(-log10(b_plasma_e107_markers$p_val_adj) > 5 & abs(b_plasma_e107_markers$avg_log2FC) >= 1, 
                                       ifelse(b_plasma_e107_markers$avg_log2FC >= 1 , 'Up' ,'Down'),' ')
writexl::write_xlsx(b_plasma_e107_markers, "OVARY5P_b_plasma_e107_treatment_DEG.xlsx" )


b_plasma_e107_markers$significant <- b_plasma_e107_markers$significant %>% as.factor()
b_plasma_e107_markers$significant <- factor(b_plasma_e107_markers$significant, levels = c("Up", "Down", " "))



ggplot(b_plasma_e107_markers, aes(avg_log2FC, -log10(p_val_adj), fill = significant, color = significant)) + 
  geom_point(size = 1.5, shape = 21) +
  theme_bw() + 
  geom_hline(yintercept=5, linetype='dashed', color='red', size=0.5) +
  geom_vline(xintercept=1, linetype = 'dashed', color='red', size = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color='red', size = 0.5) +
  scale_fill_manual(values=c("red","blue","gray67")) +
  scale_color_manual(values=c("red","blue","gray67")) +
  geom_text_repel(data = b_plasma_e107_markers[b_plasma_e107_markers[,"significant"] == "Up" | 
                                                 b_plasma_e107_markers[,"significant"] == "Down" ,],
                  aes(label = genes), 
                  size = 5,
                  color = "black") + theme(legend.text = element_text(size = 15),
                                           legend.title = element_text(size = 15),
                                           axis.text = element_text(size = 15),
                                           axis.title = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA)))) 
ggsave("b_plasma_e107_markers_volcano.png", width = 8, height = 6, dpi = 300)
ggsave("b_plasma_e107_markers_volcano.pdf", width = 8, height = 6, dpi = 600)

b_plasma_e107 <- NormalizeData(object=b_plasma_e107, normalization.method="LogNormalize", scale.factor=10000)

## Identification of highly variable features
b_plasma_e107 <- FindVariableFeatures(b_plasma_e107, selection.method = "vst", nfeatures = num.features)
top10 <- head(VariableFeatures(b_plasma_e107), 10)

## Centering the data (z-scoring)
all.genes <- rownames(b_plasma_e107)
b_plasma_e107 <- ScaleData(b_plasma_e107)

## Run PCA
library(tibble)
b_plasma_e107 <- RunPCA(b_plasma_e107, npcs=50, features = VariableFeatures(b_plasma_e107))

## Determine statistically significant principal components
#### Elbow plot
ElbowPlot(b_plasma_e107, ndims=50, reduction="pca") + theme_light()

##################################
library(Seurat)
library(bbknnR)

b_plasma_e107 <- b_plasma_e107 %>% 
  RunBBKNN(run_TSNE = FALSE, n_pcs = 8, batch_key = "orig.ident") %>% 
  FindClusters(resolution = 0.1, graph.name= "bbknn") 

celltype.id <- c("0"="Plasma cells",  "1"="B cells", "2"= "Cycling")

identical(length(celltype.id), length(unique(b_plasma_e107$seurat_clusters)))
b_plasma_e107@meta.data$celltype <- plyr::mapvalues(b_plasma_e107@meta.data$seurat_clusters, from=names(celltype.id), to=celltype.id)
b_plasma_e107$celltype %>% table

Idents(b_plasma_e107) <- "celltype"
b_plasma_e107_Markers = FindAllMarkers(object=b_plasma_e107, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

b_plasma_e107_Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
b_plasma_e107 <- ScaleData(b_plasma_e107, features = top10$gene)
av_b_plamsa <-AverageExpression(b_plasma_e107, return.seurat = TRUE)

DoHeatmap(object=av_b_plamsa, features=top10$gene, draw.lines = FALSE, angle = 45, group.colors = celltype.colv,size = 4) + theme(text = element_text(size = 12)) + NoLegend()
ggsave(file="b_plasma_e107_average_heatmap.png", width=3, height=6, dpi = 300)
ggsave(file="b_plasma_e107_average_heatmap.pdf", width=3, height=6, dpi = 600)


pt.size <- 2
p1 <- DimPlot(b_plasma_e107, reduction = "umap", group.by = "celltype", pt.size= pt.size, cols = celltype.colv) + 
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
ggsave("dimplot_celltype_b_plasma_e107.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("dimplot_celltype_b_plasma_e107.pdf", width = 3, height = 3, dpi = 600)


legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_b_plasma_e107_legend.png", width = 2, height = 2, dpi = 300)
ggsave("dimplot_celltype_b_plasma_e107_legend.pdf", width = 2, height = 2, dpi = 600)

p2 <- DimPlot(b_plasma_e107, reduction = "umap", group.by = "orig.ident", pt.size= pt.size, cols = orig.ident.colv) + 
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
ggsave("b_plasma_e107_dimplot_celltype_orig.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("b_plasma_e107_dimplot_celltype_orig.pdf", width = 3, height = 3, dpi = 600)

legend <- get_legend(p2 + theme(legend.background = element_blank(),
                                legend.title = element_blank()) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend, xlim = c(1,1) )
ggsave("b_plasma_e107_dimplot_celltype_orig_legend.png", width = 1.5, height = 3.5, dpi = 300)

b_plasma_e107$chemo <- factor(b_plasma_e107$chemo, levels = c("Pre-treatment", "Post-treatment"))
p6 <- DimPlot(b_plasma_e107, reduction = "umap", group.by = "chemo", pt.size=pt.size, cols = c("#00B0F0", "#00B050")) + 
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
ggsave("b_plasma_e107_dimplot_total_treatment.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("b_plasma_e107_dimplot_total_treatment.pdf", width = 3, height = 3, dpi = 600)

legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("b_plasma_e107+dimplot_total_treatment_legend.png", width = 3, height = 1, dpi = 300)




Input.df <- b_plasma_e107@meta.data
Input.df$celltype <- gdata::drop.levels(Input.df$celltype)

Input.df$orig.ident %>% table
Input.df$orig.ident <- factor(Input.df$orig.ident, levels = c("T70",
                                                              "T04", "T110",
                                                              "T92", "T66",
                                                              "T17", "T18",
                                                              "T39", "T34"))

Input.df$celltype <- factor(Input.df$celltype, levels = c("Plasma cells", "B cells", "Cycling"))
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

ggsave("seurat_clusters_proportion_b_plasma_e107.png", width = 6, height = 4, dpi = 300)
ggsave("seurat_clusters_proportion_b_plasma_e107.pdf", width = 6, height = 4, dpi = 600)


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
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 12),   
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size=14)) + 
  facet_wrap(~ celltype, ncol = 3) +
  scale_fill_manual(values = chemo_color) + expand_limits(y=c(0, 100)) +
  scale_color_manual(values = chemo_color) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 5,label.y = 90)

ggsave("b_plasma_e107_cluster_boxplot.png", width = 6, height = 3, dpi = 300)
ggsave("b_plasma_e107_cluster_boxplot.pdf", width = 6, height = 3, dpi = 600)


b_plasma_e107_Markers = FindAllMarkers(object=b_plasma_e107, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

b_plasma_e107_Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
b_plasma_e107 <- ScaleData(b_plasma_e107, features = top10$gene)

av_b_plasma_e107 <-AverageExpression(b_plasma_e107, return.seurat = TRUE)

DoHeatmap(object=av_b_plasma_e107, features=top10$gene, draw.lines = FALSE, size = 4) + theme(text = element_text(size = 12)) + NoLegend()

ggsave(file="b_plasma_e107_heatmap.png", width=3, height=7, dpi = 300) 


saveRDS(b_plasma_e107, "OVARY5P_final_b_plasma_e107.rds")



