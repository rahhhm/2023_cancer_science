library(ggplot2)
library(Seurat)
library(ggpubr)
library(gdata)
library(tibble)

cd4t.wd <- paste0(tnk.path, "CD4T/")
dir.create(cd4t.wd)
setwd(cd4t.wd)

Idents(tnk) <- "celltype"
tnk$celltype %>% table
cd4 <- subset(tnk, idents = c("CD4 Tn", "CD4 Tm", "CD4 T ISG", "CD4 Tfh", "CD4 Treg"))

p1 <- DimPlot(cd4, reduction = "umap", group.by = "celltype", pt.size= 2, cols = celltype.colv)  + labs(title = "") +
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
ggsave("Dimplot_celltype_cd4T.png", width = 5.5, height = 5.5, dpi = 300)
ggsave("Dimplot_celltype_cd4T.pdf", width = 5.5, height = 5.5, dpi = 600)


legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("Dimplot_celltype_cd4T_legend.png", width = 1.2, height = 2, dpi = 300)
ggsave("Dimplot_celltype_cd4T_legend.pdf", width = 1.2, height = 2, dpi = 600)



####barplot by celltype####
cd4$celltype <- gdata::drop.levels(cd4$celltype)
cd4$celltype <- factor(cd4$celltype, levels = c("CD4 Tn", "CD4 Tm", "CD4 T ISG", "CD4 Tfh",
                                                "CD4 Treg"))

Input.df <- cd4@meta.data
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
Input2$celltype <- gdata::drop.levels(Input2$celltype)

Input2$celltype <- factor(Input2$celltype, levels = c("CD4 Tn", "CD4 Tm", "CD4 T ISG", "CD4 Tfh", "CD4 pre_Treg",
                                                "CD4 Treg"))
table(Input2$celltype)
Input2 <- Input2 %>% filter(orig.ident != "T39")

count_input <- unique(Input2[,c(1,4)])
Input2$orig.ident %>% table

p1 <- ggplot(Input2, aes(x=orig.ident, y=Percent, fill=celltype)) + 
  geom_bar(stat="identity",position = position_fill(reverse = TRUE), color="black")  +
  labs(x=NULL, y="Proportion") + coord_flip() + 
  scale_fill_manual(values = celltype.colv) +
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

ggsave("celltype_proportion_cd4.png", width = 6, height = 4, dpi = 300)
ggsave("celltype_proportion_cd4.pdf", width = 6, height = 4, dpi = 600)


####barplot by clusters####

Input.df <- cd4@meta.data
Input.df$seurat_clusters <- gdata::drop.levels(Input.df$seurat_clusters)

Input.df$orig.ident %>% table
Input.df$orig.ident <- factor(Input.df$orig.ident, levels = c("T107", "T70",
                                                              "T04", "T110",
                                                              "T92", "T66",
                                                              "T17", "T18",
                                                              "T39", "T34"))
nCells2 <- dplyr::count (Input.df, orig.ident)

pre_Counted2 <- Input.df %>% group_by(orig.ident) %>% dplyr::count(seurat_clusters)
Counted2 <- pre_Counted2 %>% ungroup %>% tidyr::complete(seurat_clusters, orig.ident,  fill = list(n = 0)) # with zero-values


Input2 <- merge (x=Counted2, y=nCells2, by = "orig.ident", all.x=TRUE)
Input2$Percent <- Input2$n.x / Input2$n.y * 100

table(Input2$seurat_clusters)
Input2 <- Input2 %>% filter(orig.ident != "T39")

count_input <- unique(Input2[,c(1,4)])
Input2$orig.ident %>% table

p1 <- ggplot(Input2, aes(x=orig.ident, y=Percent, fill=seurat_clusters)) + 
  geom_bar(stat="identity",position = position_fill(reverse = TRUE), color="black")  +
  labs(x=NULL, y="Proportion") + coord_flip() +
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

ggsave("seurat_clusters_proportion_cd4.png", width = 6, height = 4, dpi = 300)

cd4$celltype <- gdata::drop.levels(cd4$celltype)
cd4$celltype <- factor(cd4$celltype, levels = c("CD4 Treg", "CD4 Tfh", "CD4 Tm", "CD4 T ISG", "CD4 Tn"))

cd4@meta.data %>% colnames
Input.df <- cd4@meta.data
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

Input4$chemo <- factor(Input4$chemo, levels =c("Pre-treatment", "Post-treatment"))            
ggplot(Input4, aes(chemo, Percent, fill = chemo)) + theme_bw() + 
  geom_boxplot(width = 0.5) +
  geom_point(aes(color = chemo), size = 1.5, shape = 19) + geom_jitter(width = 0) + NoLegend() +
  labs(y = "Proportion of CD4 T cells")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.text.x = element_text(angle = 50, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 15),   
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size=14)) + 
  facet_wrap(~ celltype, ncol = 6) +
  scale_fill_manual(values = (chemo_color)) + expand_limits(y=c(65)) +
  scale_color_manual(values = (chemo_color)) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4,label.y = 67)

ggsave(paste0("proportion_boxplot_cd4_39,107_2.png"), width = 8, height = 4, dpi = 300)
ggsave(paste0("proportion_boxplot_cd4_39,107_2.pdf"), width = 8, height = 4, dpi = 600)

# ggsave(paste0("proportion_boxplot_cd4_39_2.png"), width = 12, height = 5, dpi = 300)
# #####scatterplot####
# 
# 
# 
# markers <- list()
# markers$cytotoxic <- c("IFNG", "PRF1", "NKG7", "GZMK", "GZMA", "GNLY")
# markers$dysfunction <- c("LAG3", "PDCD1", "HAVCR2", "TIGIT", "ITGAE", "TOX", "ENTPD1")
# markers$naive <- c("CCR7", "IL7R", "TCF7", "SELL", "GPR183", "LEF1")
# markers$regulatory <- c("IL2RA","TNFRSF9", "FOXP3", "CTLA4")
# markers$ifn <- c("IFI27", "MX1", "OAS2", "ISG15", "IFIT3", "STAT1")
# 
# ######
# cd4 <- AddModuleScore(cd4, features = markers)
# 
# df <- cd4@meta.data
# df %>% colnames
# p1 <- ggplot(df, aes(x = Cluster5, y = Cluster2, color = celltype)) + 
#   geom_point(shape =19, size = 2, alpha = 1) + labs(x = "ifn", y= "ex") +
#   scale_color_manual(values = celltype.colv) + theme_bw() + 
#   theme(aspect.ratio = 1,
#         axis.title = element_text(size = 15),
#         legend.position = "bottom", 
#         legend.title = element_blank(),
#         legend.text = element_text(size = 17),
#         axis.text = element_text(size = 15))+ 
#   scale_x_continuous(limits = c(-1, 2)) + 
#   scale_y_continuous(limits = c(-1, 2)) + 
#   guides(colour = guide_legend(override.aes = list(size=4), ncol = 3))
# p1
# ggsave("cd4_score_scatterplot.png", width = 5, height = 6, dpi = 300)
# 
# p1 <- ggplot(df, aes(x = Cluster3, y = Cluster2, color = celltype)) + 
#   geom_point(shape =19, size = 2, alpha = 1) + labs(x = "IFN_score", y= "Regulatory_score") +
#   scale_color_manual(values = celltype.colv) + theme_bw() + 
#   theme(aspect.ratio = 1,
#         axis.title = element_text(size = 15),
#         legend.position = "bottom", 
#         legend.title = element_blank(),
#         legend.text = element_text(size = 17),
#         axis.text = element_text(size = 15))+ 
#   scale_x_continuous(limits = c(-1, 2)) + 
#   scale_y_continuous(limits = c(-1, 2)) + 
#   guides(colour = guide_legend(override.aes = list(size=4), ncol = 3))
# p1
# ggsave("cd4_score_scatterplot_IFN.png", width = 5, height = 6, dpi = 300)
# 
# p2 <-  ggplot(df, aes(x = Cluster1, y = Cluster2, fill = chemo)) + 
#   geom_point(shape =16, size = 2, alpha = 1 ) + labs(x = "Naive_score", y= "Regulatory_score") +
#   scale_fill_manual(values =  c("#8fffff", "#00b050")) + theme_bw() + 
#   theme(aspect.ratio = 1,
#         axis.title = element_text(size = 15), 
#         legend.position = "bottom", 
#         legend.title = element_blank(),
#         legend.text = element_text(size = 17),
#         axis.text = element_text(size = 15)) + 
#   scale_x_continuous(limits = c(-1, 2)) + 
#   scale_y_continuous(limits = c(-1, 2)) + 
#   guides(colour = guide_legend(override.aes = list(size=4)))
# 
# p1 + p2 
# ggsave("cd4_score_scatterplot.png", width =      , height = 12, dpi = 300)
# 
# library(ggplot2)
# library(hrbrthemes)
# library(dplyr)
# library(tidyr)
# library(viridis)
# 
# 
# df2 <- df %>% filter(orig.ident != "T39")
# df3 <- df2 %>% filter(orig.ident != "T107")
# df_post <- df %>%  filter(seurat_clusters == 5 | seurat_clusters == 9)
# df_naive <- df
# 
# p1 <- ggplot(data=df_naive, aes(x=Cluster1,group=seurat_clusters, fill=seurat_clusters)) +
#   geom_density(adjust=2, position="fill") + labs(x= "Naive score") +
#   theme_ipsum() 
# 
# p2 <- ggplot(data=df_post, aes(x=Cluster1,group=seurat_clusters, fill=seurat_clusters)) + 
#   geom_density(adjust=2, position="fill") + labs(x= "Naive score") +
#   theme_ipsum()
# 
# p1 + p2 +patchwork::plot_layout(ncol = 1)
# 
# ggsave("density_celltype_naive_plot.png", width =  6, height = 6, dpi = 300)
# 
# 
# p1 <- ggplot(data=df_naive, aes(x=Cluster2,group=celltype, fill=celltype)) +
#   geom_density(adjust=1.5, position="fill") + labs(x= "Regulatory score") +
#   theme_ipsum() + scale_fill_manual(values = celltype.colv)
# 
# p2 <- ggplot(data=df_post, aes(x=Cluster2,group=celltype, fill=celltype)) +
#   geom_density(adjust=1.5, position="fill") + labs(x= "Regulatory score") +
#   theme_ipsum() + scale_fill_manual(values = celltype.colv)
# 
# 
# p1 + p2 +patchwork::plot_layout(ncol = 1)
# 
# ggsave("density_celltype_regulatory_plot.png", width =  6, height = 6, dpi = 300)

####monocle3####
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(viridis)

# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3
target.obj <- cd4

# ...1 Convert to cell_data_set object ------------------------
Idents(target.obj) <- "seurat_clusters"
cds_4 <- as.cell_data_set(target.obj)
cds_4

# to get cell metadata
colData(cds_4)
# to gene metdata
fData(cds_4)
rownames(fData(cds_4))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds_4)$gene_short_name <- rownames(fData(cds_4))

# to get counts
counts(cds_4)

# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds_4@colData@rownames)))
names(reacreate.partition) <- cds_4@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds_4@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- target.obj@active.ident
cds_4@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds_4@int_colData@listData$reducedDims$UMAP <- target.obj@reductions$umap@cell.embeddings


target.obj$celltype

# ...3. Learn trajectory graph ------------------------
cds_4 <- learn_graph(cds_4, use_partition = TRUE)
root_group <- colnames(target.obj[,target.obj$seurat_clusters == "6"])
cds_4 <- order_cells(cds_4, reduction_method = 'UMAP', root_cells = root_group)
cd4$pseudotime <- pseudotime(cds_4)

p1 <- plot_cells(cds_4,
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
ggsave("monocle3_celltype_CD4T.png", width = 6, height = 5.5, dpi = 300)
ggsave("monocle3_celltype_CD4T.pdf", width = 5.5, height = 5.5, dpi = 300)


library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(viridis)

modulated_genes_cd4 <- graph_test(cds_4, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes_cd4, q_value == 0 & morans_I > 0.23))
genes <- genes[-grep(pattern = "RP", genes)]

pt.matrix <- normalized_counts(cds_4, norm_method = "log")[match(genes,rownames(rowData(cds_4))),order(pseudotime(cds_4))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds_4, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
colnames(pt.matrix) <- cd4$pseudotime[normalized_counts(cds_4, norm_method = "log")[match(genes,rownames(rowData(cds_4))),order(pseudotime(cds_4))] %>% colnames()]

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
col_fun = colorRamp2(c(colnames(pt.matrix) %>% as.numeric()), c(viridis::plasma(3272)))
ha = HeatmapAnnotation(pseudotime = colnames(pt.matrix) %>% as.numeric(), col = list(pseudotime = col_fun))

hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 10),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = ha)

png("CD8_monocle_heatmap.png",width=6,height=8,units="in",res=1200)
print(hthc)
dev.off()

















