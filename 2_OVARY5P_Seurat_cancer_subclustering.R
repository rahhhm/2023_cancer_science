library(Seurat)
library(dplyr)
library(bbknnR)
library(stringr)

sample.list <- list.files("d:/project_directory_name/ovary_project/OVARY_5P/", pattern = "5P")

cancer$chemo <- cancer$orig.ident
new.meta.ids <- c('Post-treatment','Post-treatment','Post-treatment',
                  'Pre-treatment','Pre-treatment','Pre-treatment','Pre-treatment',
                  'Pre-treatment','Post-treatment','Pre-treatment')
cancer$chemo <- plyr::mapvalues(cancer$chemo, from = sample.list, to = new.meta.ids )
table(cancer$orig.ident, cancer$chemo)


cancer$orig.ident <- str_split_fixed(cancer$orig.ident, "-", 3)[,2]
cancer$orig.ident <- factor(cancer$orig.ident, levels = rev(c("T107", "T70",
                                                            "T04", "T110",
                                                            "T92", "T66",
                                                            "T17", "T18",
                                                            "T39", "T34")))

####cancer####
label <- "epi"
epi.path <- paste0(wd.path,"epithelial cells/")
setwd(epi.path)

Idents(ovary) <- "celltype"
cancer <- subset(ovary, idents = "Epithelial cells")

num.features <- 2000

cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = num.features)
top10 <- head(VariableFeatures(cancer), 10)

## Centering the data (z-scoring)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cancer <- CellCycleScoring(cancer, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cancer <- ScaleData(cancer, vars.to.regress = c("s.genes", "G2M.Score"), features = VariableFeatures(cancer))

#, vars.to.regress = "percent.mt"

## Run PCA
library(tibble)
cancer <- RunPCA(cancer, npcs=50)

## Determine statistically significant principal components
#### Elbow plot
ElbowPlot(cancer, ndims=50, reduction="pca")

# cancer <- cancer %>% 
#   RunBBKNN(run_TSNE = FALSE, n_pcs = 10, batch_key = "orig.ident") %>% 
#   FindClusters(resolution = 1.0, graph.name= "bbknn")

cancer <- FindNeighbors(cancer, dims = 1:20) %>% 
  FindClusters(resolution = 1.5) %>% 
  RunUMAP(dims = 1:20)

#####Dimplot of cancer cells####
p2 <- DimPlot(cancer, reduction = "umap", group.by = "orig.ident", pt.size= pt.size, cols = orig.ident.colv) + 
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
ggsave(paste0(label,"_dimplot_celltype_orig.png"), width = 6, height = 6, dpi = 300)
ggsave(paste0(label,"_dimplot_celltype_orig.pdf"), width = 6, height = 6, dpi = 600)

legend <- get_legend(p2 + theme(legend.background = element_blank(),
                                legend.title = element_blank()) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend, xlim = c(1,1) )
ggsave(paste0(label,"_dimplot_celltype_orig_legend.png"), width = 2, height = 6, dpi = 300)

p6 <- DimPlot(cancer, reduction = "umap", group.by = "chemo", pt.size=pt.size, cols = c("#00B0F0", "#00B050")) + 
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
AugmentPlot(p6)
ggsave(paste0(label,"dimplot_total_tratment.png"), width = 6, height = 6, dpi = 300)
ggsave(paste0(label,"dimplot_total_tratment.pdf"), width = 6, height = 6, dpi = 600)

legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave(paste0(label,"dimplot_total_tratment_legend.png"), width = 4, height = 1, dpi = 300)
ggsave(paste0(label,"dimplot_total_tratment_legend.pdf"), width = 4, height = 1, dpi = 600)




library(GSVA)
library(msigdbr)
library(msigdb)
library(pheatmap)

collections <- msigdbr_collections()
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", category = "H")

hallmarks_list<- hallmark_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)


hallmarks_name<- names(hallmarks_list)
for (i in 1 : length(hallmarks_list)) {
  hallmarks_list[[i]] <- unique(hallmarks_list[[i]])
}
names(hallmarks_list) <- hallmarks_name

cancer <- ScaleData(cancer, features = c(unlist(hallmarks_list) %>% as.character() %>% unique()))
avg <- AverageExpression(cancer, group.by = "orig.ident", slot = 'scale.data')
avg <- avg$RNA
avg <- avg[,-10]
colnames(avg) 

avg <- avg[,rev(c(8,1,2,9,7,3,4,6,5))]
avg <- avg[,rev(c(9,1,3,2,10,8,4,5,7,6))]



#check gene list in matrix
for (i in 1:length(hallmarks_list)) {
  count = hallmarks_list[[i]] %in% rownames(avg) %>% table
  print(count)
}

gsva_results <- GSVA::gsva( 
  avg %>% as.matrix(), 
  hallmarks_list, 
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)
rownames(gsva_results) <- gsub("HALLMARK_", "", row.names(gsva_results))

my_sample_col1 <- data.frame(Chemotherapy = rep(c("Pre-treatment", "Post-treatment"), c(6,4)))
row.names(my_sample_col1) <- colnames(avg)

group_row <- data.frame(group_name = rep(c("Group_1", "Group_2", "Group_3", "Group_4"),c(7, 22,7,14)))
row.names(group_row) <- c(rownames(gsva_results)[c(31,13,18,41,46,32,33,
                                                   35,28,21,37,9,22,19,3,48,50,1,17,15,16,8,6,39,49,12,40,30,43,
                                                   4,14,47,44,20,5,34)],
                          
                          rownames(gsva_results)[-c(31,13,18,41,46,32,33,
                                                    35,28,21,37,9,22,19,3,48,50,1,17,15,16,8,6,39,49,12,40,30,43,
                                                    4,14,47,44,20,5,34)])


group.color <- colorRampPalette(colors = brewer.pal(n = 8, name = "Spectral"))(4)

my_color <- list(
  Chemotherapy = c("Pre-treatment" = "#00B0F0", 
                   "Post-treatment" = "#00B050"),
  group_name = c("Group_1" = group.color[1], 
                 "Group_2" = group.color[2], 
                 "Group_3" = group.color[3], 
                 "Group_4" = group.color[4]))

pheatmap(gsva_results, legend_breaks = , 
         annotation_row = group_row,
         annotation_col = my_sample_col1,
         annotation_colors = my_color, 
         cutree_rows = 4,
         cluster_cols = FALSE, 
         cluster_rows = TRUE, 
         scale = 'none', filename = "tumor_gsva_heatmap.pdf",
         annotation_legend = TRUE, width = 8, height = 10,
         )
dev.off()

library(VISION)
library(ggpubr)
library(scCustomize)

signatures <- c("h.all.v2022.1.Hs.symbols.gmt")
cancer2 <- cancer
cancer@meta.data %>% colnames()
cancer2@meta.data <- cancer@meta.data[,c(42,1)]
vis <- Vision(cancer2,
               signatures = signatures)

# Set the number of threads when running parallel computations
# On Windows, this must either be omitted or set to 1
options(mc.cores = 1)

vis <- calcSignatureScores(vis)
rownames(vis@SigScores) <- vis@exprData %>%  colnames()
vis <- clusterSigScores(vis)
vis <- analyzeLocalCorrelations(vis)
vis
viewResults(vis)
head(getSignatureAutocorrelation(vis))

summary(cancer@meta.data$HALLMARK_WNT_BETA_CATENIN_SIGNALING)

saveRDS(vis, "vision_hallmark_result.rds")

vis <- readRDS("vision_hallmark_result.rds")
viewResults(vis)

sig <- vis@SigScores %>%  as.data.frame()
cancer@meta.data$EPITHELIAL_MESENCHYMAL_TRANSITION <- sig$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION 

p4 <- FeaturePlot(cancer, features = "EPITHELIAL_MESENCHYMAL_TRANSITION", pt.size = pt.size, max.cutoff = 1.3) + 
  scale_colour_gradientn(colours = c("white","#ebd4cb", "#da9f93", "#b6465f", "#890620", "#2c0703", "black")) +
  theme(aspect.ratio = 1,
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

AugmentPlot(p4)
ggsave("featureplot_cancer_emt.png", width = 6, height = 6, dpi = 300)
ggsave("featureplot_cancer_emt.pdf", width = 6, height = 6, dpi = 600)

legend <- get_legend(p4)

ggdraw(legend, xlim = c(1,1))
ggsave("featureplot_cancer_emt_legend.png", width = 2, height = 2, dpi = 300)
ggsave("featureplot_cancer_emt_legend.pdf", width = 2, height = 2, dpi = 600)

cancer@meta.data$E2F_TARGETS <- sig$HALLMARK_E2F_TARGETS 

p5 <-  FeaturePlot(cancer, features = "E2F_TARGETS", pt.size = pt.size) + 
  scale_colour_gradientn(colours = c("white","#ebd4cb", "#da9f93", "#b6465f", "#890620", "#2c0703", "black")) +
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
ggsave("featureplot_cancer_e2f.png", width = 6, height = 6, dpi = 300)
ggsave("featureplot_cancer_e2f.pdf", width = 6, height = 6, dpi = 600)

legend <- get_legend(p5)

ggdraw(legend, xlim = c(1,1))
ggsave("featureplot_cancer_e2f_legend.png", width = 2, height = 2, dpi = 300)
ggsave("featureplot_cancer_e2f_legend.pdf", width = 2, height = 2, dpi = 600)

cellcycle_color <- c("G1" = "#00a352", "S" = "#1c5d99", "G2M" = "#da1b2b")
####cell cycle ####
p6 <- DimPlot(cancer, reduction = "umap", group.by = "Phase", pt.size=pt.size, cols = cellcycle_color) + 
  labs(title = "Cell cycle phase") +
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
ggsave("Dimplot_cell_cycle_phase.png", width = 6, height = 6, dpi = 300)
ggsave("Dimplot_cell_cycle_phase.pdf", width = 6, height = 6, dpi = 600)

legend <- get_legend(p6)

ggdraw(legend, xlim = c(1,1))
ggsave("Dimplot_cell_cycle_phase_legend.png", width = 2, height = 2, dpi = 300)
ggsave("Dimplot_cell_cycle_phase_legend.pdf", width = 2, height = 2, dpi = 600)

library(ggpubr)
cancer.meta <- cancer@meta.data
cancer.meta <- cancer.meta %>% filter(doublet_result == "Singlet")
cancer.meta %>% colnames
orig.count <- dplyr::count(cancer.meta, orig.ident)
pre_phase.count <- cancer.meta %>% group_by(orig.ident) %>% dplyr::count(Phase)
phase.count <- pre_phase.count %>% ungroup %>% tidyr::complete(Phase, orig.ident,  fill = list(n = 0)) # with zero-values
cycle_table <- merge(phase.count, orig.count, by = "orig.ident")

cycle_table$prop <- cycle_table$n.x / cycle_table$n.y * 100
colnames(cycle_table) <- c("orig.ident", "phase", "count", "total", "prop")

cycle_table <- merge(cycle_table, unique(cancer.meta[,c(1,42)]), by = "orig.ident")

cycle_table <- cycle_table %>%  filter(orig.ident != "T107")

cycle_table$chemo <- factor(cycle_table$chemo, levels = c("Pre-treatment", "Post-treatment"))

cycle_table_g1 <- cycle_table %>%  filter(phase == "G1")
cycle_table_S <- cycle_table %>%  filter(phase == "S" )
cycle_table_G2M <- cycle_table %>%  filter(phase == "G2M")

my_comparisons <- list(c("Pre-treatment", "Post-treatment"))
chemo_color <- c("#00b0f0", "#00b050")


g1 <- ggplot(cycle_table_g1, aes(x = chemo, y = prop, fill = chemo)) + theme_bw() +
  geom_boxplot()  + labs(y= "Proportion",title = "G1 phase") + scale_fill_manual(values = chemo_color) + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", size = 7, label.y = 85) + geom_jitter(width = 0) +
  theme(aspect.ratio= 1.8,
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 15),
        title = element_text(size = 17),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank()) + NoLegend() + expand_limits(y=c(25,95))

s <- ggplot(cycle_table_S, aes(x = chemo, y = prop, fill = chemo)) + theme_bw() +
  geom_boxplot()  + labs(title = "S phase") + scale_fill_manual(values = chemo_color) + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", size = 7, label.y = 55) + geom_jitter(width = 0) + 
  theme(aspect.ratio= 1.8,
        axis.title.x.bottom = element_blank(),
        title = element_text(size = 17),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank()) + NoLegend() + expand_limits(y=c(5,60))

g2m <- ggplot(cycle_table_G2M, aes(x = chemo, y = prop, fill = chemo)) + theme_bw() +
  geom_boxplot()  + labs(title = "G2M phase") + scale_fill_manual(values = chemo_color) + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", size = 7, label.y = 46) + geom_jitter(width = 0) + NoLegend() +
  theme(aspect.ratio= 1.8,
        axis.title.x.bottom = element_blank(),
        title = element_text(size = 17),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_blank()) + expand_limits(y=c(5,50))

g1 + s + g2m

ggsave("cellcycle_boxplot.png", width = 6.5, height = 4, dpi = 300)
ggsave("cellcycle_boxplot.pdf", width = 6, height = 5, dpi = 600)

cycle_table_pair <- cycle_table %>% filter(orig.ident == "T92" | orig.ident == "T110")


ggplot(cycle_table_pair, aes(x=orig.ident, y=prop, fill=phase)) + 
  geom_bar(stat="identity",position = position_fill(reverse = TRUE), color="black")  +
  labs(x=NULL, y="% Proportion") + coord_flip() + scale_fill_manual(values = cellcycle_color) +
  theme_classic() + theme(legend.position="bottom",
                          legend.direction = "horizontal",
                          legend.key.size =  unit(1,"cm"),
                          legend.title = element_blank(),
                          legend.text = element_text(colour="black", size=12),
                          legend.background = element_blank(),
                          axis.text.y = element_text(size = 20),
                          axis.text.x = element_text(size = 15), 
                          axis.title.x = element_text(size = 15)) 

ggsave("cancer_pair_celltype_barplot.png", width = 4, height = 3, dpi = 300)
ggsave("cancer_pair_celltype_barplot.pdf", width = 5, height = 3, dpi = 600)

#########DEG analysis_GSEA##########


library(dplyr)
library(ggrepel)
library(fgsea)
library(reshape2)
library(ggrepel)
library(ggtext)
library(glue)

Idents(cancer) <- "chemo"
clustering <- cancer@active.ident
cellsIn <- names(clustering[clustering == "Post-treatment"])
cellsOut <- names(clustering[clustering == "Pre-treatment"])
cancer<- ScaleData(cancer)

cancer_markers <- FindMarkers(cancer, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)
rownames(cancer_markers) <- cancer_markers$genes
### MsigDB doesn't play well with mouse genes, using msigdbr for help
ranks <- cancer_markers$avg_log2FC
names(ranks) <- rownames(cancer_markers)
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", category = "H")

msigdbr_list = split (x=hallmark_gene_sets$gene_symbol, f=hallmark_gene_sets$gs_name)

fgseaRes <- fgsea(pathways=msigdbr_list, stats=ranks, minSize=5, maxSize=1500)
sig_fgseaRes <- fgseaRes %>% filter(padj < 0.05)
sig_fgseaRes$pathway <- gsub("HALLMARK_", "", sig_fgseaRes$pathway)

sig_fgseaRes$chemo <- ifelse(sig_fgseaRes$NES > 0, "Post-treatment", "Pre-treatment")
chemo_color <- c("Pre-treatment" = "#00B0F0","Post-treatment" = "#00B050")

bold.labels <- ifelse(arrange(sig_fgseaRes, desc(NES))$pathway %in% signif_list, yes = "bold", no = "plain") %>%  rev()
size.labels <- ifelse(arrange(sig_fgseaRes, desc(NES))$pathway %in% signif_list, yes = 12, no = 10) %>%  rev()

ggplot(sig_fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= chemo)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal(base_size = 15) + scale_fill_manual(values = chemo_color) + 
  theme(axis.text.y = element_text(face = bold.labels, size = size.labels))


ggsave("fgsea_result.png", width = 9, height = 5, dpi = 300)
ggsave("fgsea_result.pdf", width = 7, height = 5, dpi = 600)




cancer_markers$genes <- rownames(cancer_markers)
cancer_markers$significant <- ifelse(cancer_markers$p_val_adj < 0.05 & abs(cancer_markers$avg_log2FC) >= 1, 
                                     ifelse(cancer_markers$avg_log2FC >= 1 , 'Up' ,'Down'),' ')


cancer_markers$significant <- cancer_markers$significant %>% as.factor()
cancer_markers$significant <- factor(cancer_markers$significant, levels = c("Up", "Down", " "))

unique(melt(sig_fgseaRes$leadingEdge))$value[unique(melt(sig_fgseaRes$leadingEdge))$value %in% 
                                               rownames(cancer_markers %>% filter(abs(avg_log2FC) >1))]


label <- unique(melt(sig_fgseaRes$leadingEdge)$value)
cancer_markers$significant %>% table

ggplot(cancer_markers, aes(avg_log2FC, -log10(p_val_adj), fill = significant, color = significant)) + 
  geom_point(size = 1.5, shape = 21) +
  theme_bw() + 
  geom_hline(yintercept=1.3, linetype='dashed', color='red', size=0.5) +
  geom_vline(xintercept=1, linetype = 'dashed', color='red', size = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color='red', size = 0.5) +
  geom_text_repel(data = cancer_markers[label[label %in% rownames(cancer_markers %>% filter(abs(avg_log2FC) >1))],],
                  aes(label = genes), 
                  size = 4,
                  color = "black",
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))+
  scale_fill_manual(values=c("red","blue", "gray67")) +
  scale_color_manual(values=c("red","blue", "gray67")) + theme(legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA)))) 
  

ggsave("volcano_total.png", width = 8, height = 5, dpi = 300)
ggsave("volcano_total.pdf", width = 8, height = 5, dpi = 600)


write.csv(cancer_markers, "total_treatmentdeg.csv")

#####cancer_pair volcano plot####
Idents(cancer) <- "orig.ident"
cancer_pair <- subset(cancer, idents = c("T110", "T92"))

Idents(cancer_pair) <- "chemo"
clustering <- cancer_pair@active.ident
cellsIn <- names(clustering[clustering == "Post-treatment"])
cellsOut <- names(clustering[clustering == "Pre-treatment"])

cancer_pair<- ScaleData(cancer_pair)
cancer_markers2 <- FindMarkers(cancer_pair, ident.1 = cellsIn, ident.2 = cellsOut, logfc.threshold = 0, test.use = 'MAST')


### MsigDB doesn't play well with mouse genes, using msigdbr for help
ranks <- cancer_markers2$avg_log2FC
names(ranks) <- rownames(cancer_markers2)
fgseaRes2 <- fgsea(pathways=msigdbr_list, stats=ranks, minSize=5, maxSize = 1500 )
sig_fgseaRes2 <- fgseaRes2 %>% filter(padj < 0.01)

sig_fgseaRes2$pathway <- gsub("HALLMARK_", "", sig_fgseaRes2$pathway)

sig_fgseaRes2$chemo <- ifelse(sig_fgseaRes2$NES > 0, "Post-treatment", "Pre-treatment")


sig_fgseaRes2$signif <- ifelse(sig_fgseaRes2$pathway %in% rownames(filter(comparae_fgsea, signif2 == "signif")), "signif", "none")

bold.labels <- ifelse(arrange(sig_fgseaRes2, desc(NES))$pathway %in% signif_list, yes = "bold", no = "plain") %>%  rev()
size.labels <- ifelse(arrange(sig_fgseaRes2, desc(NES))$pathway %in% signif_list, yes = 12, no = 10) %>%  rev()

ggplot(sig_fgseaRes2, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill= chemo)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal(base_size = 15) + scale_fill_manual(values = chemo_color) +
  theme(axis.text.y = element_text(face = bold.labels))

ggsave("fgsea_pair_result.png", width = 9, height = 6, dpi = 300)
ggsave("fgsea_pair_result.pdf", width = 9, height = 5.5, dpi = 600)



comparae_fgsea %>% filter(signif2 == "signif") %>% rownames

cancer_markers2 <- cancer_markers2 %>% arrange(p_val_adj)
cancer_markers2$genes <- rownames(cancer_markers2)
cancer_markers2 <- cancer_markers2 %>% filter(p_val_adj < 0.05)
cancer_markers2$significant <- ifelse(abs(cancer_markers2$avg_log2FC) >= 1, 
                                      ifelse(cancer_markers2$avg_log2FC >= 1 , 'Up' ,'Down'),' ')

label2 <- unique(melt(sig_fgseaRes2$leadingEdge)$value)
cancer_markers2$significant %>% table

cancer_markers2$significant <- cancer_markers2$significant %>% as.factor()
cancer_markers2$significant <- factor(cancer_markers2$significant, levels = c("Up", "Down", " "))



ggplot(cancer_markers2, aes(avg_log2FC, -log10(p_val_adj), fill = significant, color = significant)) + 
  geom_point(size = 1.5, shape = 21) + 
  theme_bw() + 
  geom_hline(yintercept=1.3, linetype='dashed', color='red', size=0.5) +
  geom_vline(xintercept=1, linetype = 'dashed', color='red', size = 0.5) +
  geom_vline(xintercept=-1, linetype = 'dashed', color='red', size = 0.5) +
  geom_text_repel(data = cancer_markers2[label[label %in% rownames(cancer_markers2 %>% filter(abs(avg_log2FC) >1))],],
                  aes(label = genes), 
                  size = 4,
                  color = "black",
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))+
  scale_fill_manual(values=c("red","blue", "gray67")) +
  scale_color_manual(values=c("red","blue", "gray67"))+ theme(legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA)))) 

ggsave("volcano_pair_total.png", width = 8, height = 5, dpi = 300)
ggsave("volcano_pair_total.pdf", width = 8, height = 5, dpi = 600)


writexl::write_xlsx(cancer_markers2, "../../epithelial cells/pair_treatmentdeg.xlsx")


#fgseaRes, fgseaRes2
fgseaRes <- fgseaRes %>% arrange()
fgseaRes2 <- fgseaRes2 %>% arrange()


comparae_fgsea <- merge(x = fgseaRes[,c(1,6)] , y = fgseaRes2[,c(1,6)], by = "pathway")
comparae_fgsea <- column_to_rownames(comparae_fgsea, var = "pathway")
#comparae_fgsea %>%  plyr::sort(decreasing = T)
colnames(comparae_fgsea) <- c("Total", "Pair")
rownames(comparae_fgsea) <- gsub("HALLMARK_", "", rownames(comparae_fgsea))
comparae_fgsea$signif <- ifelse(abs(comparae_fgsea$Total) > 1.5 | abs(comparae_fgsea$Pair) > 1.5,
                                ifelse(abs(comparae_fgsea$Total) > 1.5 & abs(comparae_fgsea$Pair) > 1.5, "Both_signif",
                                       ifelse(abs(comparae_fgsea$Pair) > 1.5 & abs(comparae_fgsea$Total) < 1.5,"Pair_signif", "Total_signif") ), "none")

comparae_fgsea$signif2 <- ifelse(abs(comparae_fgsea$Total) > 1.5 & abs(comparae_fgsea$Pair) > 1.5, "signif", "none")

comparae_fgsea$signif2 <- factor(comparae_fgsea$signif2, levels = c("signif", "none"))
ggplot(comparae_fgsea, aes(x = Total, y = Pair)) + geom_point(aes(color = signif2)) + theme_bw() + labs(x = "NES(total)", y = "NES(pair)", title = "HALLMARK_SIGNATURES") + guides(color = guide_legend(title="padj < 0.05"), ) + theme(text = element_text(size = 17)) +
  stat_cor(method = "pearson", data = comparae_fgsea, label.y = 0.5,label.x = -2, size = 5 ) + scale_color_manual(values = c("red","black")) +  geom_smooth(method='lm', formula = 'y ~ x') 
ggsave("scatter_fgsea.png", width = 7, height = 5, dpi = 300)
ggsave("epithelial cells/scatter_fgsea.pdf", width = 7, height = 5, dpi = 600)

signif_list <- comparae_fgsea %>%  filter(signif == "Both_signif") %>%  rownames
signif_list <- comparae_fgsea %>%  filter(signif == "Both_signif") %>%  rownames
