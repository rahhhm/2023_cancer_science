library(cowplot)
library(patchwork)

cd8t.wd <- paste0(tnk.path, "CD8T/")
dir.create(cd8t.wd)
setwd(cd8t.wd)

Idents(tnk) <- "celltype"
DimPlot(tnk, label = TRUE, label.size = 7)
tnk$celltype %>% table
cd8 <- subset(tnk, idents = c("CD8 Tn/m", "CD8 Tem",
                              "CD8 T ISG", "CD8 Tex"))

pt.size =2
Idents(cd8) <- "seurat_clusters"

p1 <- DimPlot(cd8, reduction = "umap", group.by = "celltype", pt.size= 2, cols = celltype.colv)  + labs(title = "") +
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
  guides(colour = guide_legend(override.aes = list(size=7), ncol = 1, byrow = FALSE))
AugmentPlot(p1)
ggsave("Dimplot_celltype_CD8T.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("Dimplot_celltype_CD8T.pdf", width = 5.5, height = 5.5, dpi = 600)

legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("Dimplot_celltype_CD8T_legend.png", width = 1.2, height = 2, dpi = 300)
ggsave("Dimplot_celltype_CD8T_legend.pdf", width = 1.2, height = 2, dpi = 600)


## Percent bar graph (by cell.type - orig.ident)
cd8$celltype <- gdata::drop.levels(cd8$celltype)

cd8$celltype <- gdata::drop.levels(cd8$celltype)
cd8$celltype <- factor(cd8$celltype, levels = c("CD8 Tn/m", "CD8 Tem",
                                                          "CD8 T ISG", "CD8 Tex"))
Input.df <- cd8@meta.data
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
Input2 <- Input2 %>% filter(orig.ident != "T39")
count_input <- unique(Input2[,c(1,4)])


p1 <- ggplot(Input2, aes(x=orig.ident, y=Percent, fill=celltype)) + 
  geom_bar(stat="identity",position = position_fill(reverse = TRUE), color="black")  +
  labs(x=NULL, y="Proportion") + coord_flip()  + scale_fill_manual(values = celltype.colv) +
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

ggsave("celltype_proportion_cd8.png", width = 6, height = 4, dpi = 300)
ggsave("celltype_proportion_cd8.pdf", width = 6, height = 4, dpi = 600)


cd8$celltype <- factor(cd8$celltype, levels = c("CD8 Tem","CD8 Tex", "CD8 T ISG", "CD8 Tn/m"))
       
Input.df <- cd8@meta.data
Input.df %>% colnames
nCells2 <- dplyr::count (Input.df, orig.ident)
nCells2 <- merge(nCells2, Input.df[,c(1,41)] %>% unique() %>% arrange(desc(chemo)), by = "orig.ident", all.x=TRUE)
pre_Counted2 <- Input.df %>% group_by(orig.ident) %>% dplyr::count(celltype)
Counted2 <- pre_Counted2 %>% ungroup %>% tidyr::complete(celltype, orig.ident,  fill = list(n = 0)) # with zero-values
Input2 <- merge (x=Counted2, y=nCells2, by = "orig.ident", all.x=TRUE)
Input2$Percent <- Input2$n.x / Input2$n.y * 100
Input2$chemo %>% table

Input3 <- Input2 %>% filter(orig.ident != "T39")
Input4 <- Input3 %>% filter(orig.ident != "T107")

my_comparisons <- list(c("Pre-treatment", "Post-treatment"))

Input4$chemo <- factor(Input4$chemo, levels = c("Pre-treatment", "Post-treatment"))       
p1 <- ggplot(Input4, aes(chemo, Percent, fill = chemo)) + theme_bw() + NoLegend() + 
  geom_boxplot(width = 0.5) +
  geom_point(size = 1.5, shape = 19)  + 
  labs(y = "Proportion of CD8 T cells")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.text.x = element_text(angle = 50, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 15),   
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size=14)) + 
  facet_wrap(~ celltype, ncol = 4) +
  scale_fill_manual(values = chemo_color) + expand_limits(y=c(0, 90)) +
  scale_color_manual(values = chemo_color) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4,label.y = 70)


nCells2 <- dplyr::count (Input.df, orig.ident)

Input.df$celltype2 <- Input.df$celltype

celltype.id <- c("CD8 Tex" = "CD8 Tex + ISG", "CD8 Tem" = "CD8 Tem",
                 "CD8 T ISG" = "CD8 Tex + ISG", "CD8 Tn/m" = "CD8 Tn/m")

identical(length(celltype.id), length(unique(Input.df$celltype2)))
Input.df$celltype2 <- plyr::mapvalues(Input.df$celltype2, from=names(celltype.id), to=celltype.id)

nCells2 <- dplyr::count (Input.df, orig.ident)
nCells2 <- merge(nCells2, Input.df[,c(1,41)] %>% unique() %>% arrange(desc(chemo)), by = "orig.ident", all.x=TRUE)
pre_Counted2 <- Input.df %>% group_by(orig.ident) %>% dplyr::count(celltype2)
Counted2 <- pre_Counted2 %>% ungroup %>% tidyr::complete(celltype2, orig.ident,  fill = list(n = 0)) # with zero-values
Input2 <- merge (x=Counted2, y=nCells2, by = "orig.ident", all.x=TRUE)
Input2$Percent <- Input2$n.x / Input2$n.y * 100

Input3 <- Input2 %>% filter(orig.ident != "T39")
Input4 <- Input3 %>% filter(orig.ident != "T107")

my_comparisons <- list(c("Pre-treatment", "Post-treatment"))
Input5 <- Input4 %>% filter(celltype2 != "CD8 Tem" & celltype2 != "CD8 Tn/m")
Input5$chemo <- factor(Input5$chemo, levels = c("Pre-treatment", "Post-treatment"))       
p1 + ggplot(Input5, aes(chemo, Percent, fill = chemo)) + theme_bw() + NoLegend() + 
  geom_boxplot(width = 0.5) +
  geom_point(, size = 1.5, shape = 19) +
  labs(y = "")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),   
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size=14)) + 
  facet_wrap(~ celltype2, ncol = 4) +
  scale_fill_manual(values = chemo_color) + expand_limits(y=c(0, 90)) +
  scale_color_manual(values = chemo_color) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4,label.y = 75)  + 
  plot_layout(widths = c(4, 1))

ggsave(paste0("proportion_boxplot_cd8_39,107_2.pdf"), width = 8, height = 4, dpi = 600)

# cd8 <- NormalizeData(object=cd8, normalization.method="LogNormalize", scale.factor=10000)
# 
# ## Identification of highly variable features
# cd8 <- FindVariableFeatures(cd8, selection.method = "vst", nfeatures = num.features)
# top10 <- head(VariableFeatures(cd8), 10)
# 
# ## Centering the data (z-scoring)
# all.genes <- rownames(cd8)
# cd8 <- ScaleData(cd8, features = all.genes)
# 
# cd8 <- RunPCA(cd8, npcs=50, features = VariableFeatures(cd8))
# 
# ElbowPlot(cd8, ndims=50, reduction="pca") + theme_light()
# 
# cd8 <- cd8 %>% 
#   RunBBKNN(run_TSNE = FALSE, n_pcs = 15, batch_key = "orig.ident") %>% 
#   FindClusters(resolution = 1.5, graph.name= "bbknn")
# 
# DimPlot(cd8, label = TRUE, label.size = 7)


####monocle3####
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(viridis)
# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3
target.obj <- cd8

# ...1 Convert to cell_data_set object ------------------------
Idents(target.obj) <- "seurat_clusters"
cds <- as.cell_data_set(target.obj)
cds

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)

# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- target.obj@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- target.obj@reductions$umap@cell.embeddings

# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = TRUE)
# ...4. Order the cells in pseudotime -------------------
root_group <- colnames(target.obj[,target.obj$seurat_clusters == "12"])
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = root_group )
cd8$pseudotime2 <- pseudotime(cds)

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime", 
                 label_cell_groups=FALSE, show_trajectory_graph = TRUE, scale_to_range = F,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 label_branch_points=FALSE,
                 graph_label_size= 5,
                 rasterize = TRUE,
                 group_label_size = 5, cell_size = 0.9,
                 trajectory_graph_segment_size = 2, trajectory_graph_color = "black") +
  theme(aspect.ratio= 1, 
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        rect = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        legend.text = element_text(size =15),
        legend.title = element_text(size =15),
        axis.title = element_blank()) + scale_color_gradientn(colors=viridis::plasma(10)[2:10]) +
  labs(colour= "Pseudo-time")
  

p1 
ggsave("monocle3_celltype_CD8T.png", width = 6, height = 5.5, dpi = 300)
ggsave("monocle3_celltype_CD8T.pdf", width = 5.5, height = 5.5, dpi = 300)

legend <- get_legend(p1) 

ggdraw(legend)
ggsave("monocle3_celltype_CD8T_legend.png", width = 1.2, height = 2, dpi = 300)

#1000 500


FeaturePlot(cd8, "pseudotime2", pt.size = 0.1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(aspect.ratio= 1,
        legend.title = element_text(size = 16),
        legend.text = element_text(size=15),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13), 
        title = element_text(size = 15))

# cells ordered by monocle3 pseudotime
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(viridis)
#subset(modulated_genes, q_value == 0 & morans_I > 0.2)
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(modulated_genes)
genes <-  c("LMNA", "IL7R", "CCR7", "GZMK", 
            "TNFRSF9", "CTLA4", "KRT86", "PHLDA1", 
            "LAYN", "PLPP1", "CXCL13", "ENTPD1", 
            "HOPX", "HAVCR2", "TIGIT",
           "RGS1", "GZMB", "CD63", "TRAV21", 
           "LAG3", "HLA-DRB1", "HLA-DRA", "NKG7", "CCL5", 
           "GZMA", "OAS1", "IFIT1", "IFIT3", "RSAD2", "ISG15", "MX1") %>%  as.character()

pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
colnames(pt.matrix) <- cd8$pseudotime2[normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))] %>% colnames()]

# #K means with 6 groups
# htkm <- Heatmap(
#   pt.matrix,
#   name                         = "z-score",
#   col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
#   show_row_names               = TRUE,
#   show_column_names            = FALSE,
#   row_names_gp                 = gpar(fontsize = 6),
#   km = 6,
#   row_title_rot                = 0,
#   cluster_rows                 = TRUE,
#   cluster_row_slices           = FALSE,
#   cluster_columns              = FALSE)
#print(htkm)

#Ward.D2 Hierarchical Clustering
ht_opt$message = FALSE
col_fun = colorRamp2(c(colnames(pt.matrix) %>% as.numeric()), c(viridis::plasma(2860)))
ha = HeatmapAnnotation(pseudotime = colnames(pt.matrix) %>% as.numeric(), 
                       col = list(pseudotime = col_fun), 
                       annotation_legend_param = list(
                         pseudotime = list(#legend_direction = "horizontal", 
                                           legend_height = unit(4, "cm"),ncol = 1)))

pdf("CD8_monocle_heatmap.pdf",width=6,height=6)

hthc <- Heatmap(
  pt.matrix, 
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 13),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = ha,
  heatmap_legend_param = list(ncol = 1, legend_height = unit(6, "cm")) )


draw(hthc, heatmap_legend_side="right", annotation_legend_side="right",
      merge_legend = TRUE)
dev.off()




