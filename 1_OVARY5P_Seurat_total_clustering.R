library(SoupX)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)  
library(DoubletFinder)
library(bbknnR)
library(reshape2)
library(RColorBrewer)
library(ggrastr)


input.path <- "D:/project_directory_name/ovary_project/OVARY_5P/"
label <- "OVARY5P"
mother.path <- "D:/project_directory_name/ovary_project/Seurat/"
wd.path <- paste0(mother.path, "final/")
ifelse(!dir.exists(wd.path), dir.create(wd.path),"directory already exists")
setwd(wd.path)


# for variable genes
num.features = 3000

# Cell QC cutoff
nFeature.min = 400
nFeature.max = 9000
percent.mito.max = 20
nCount_RNA.min = 600
nCount_RNA.max = 100000
doublet_table_total <-c()
min.cells = 3



for (sam in sample.list){
  sample.path <- paste0(wd.path, "soupx_doublet_result/")
  ifelse(!dir.exists(sample.path),dir.create(sample.path),"directory already exists")
  sc = load10X(paste0(input.path, sam, "/outs/count/"))
  sc = autoEstCont(sc)
  out = adjustCounts(sc)
  colnames(out) <- paste0(sam, "_", colnames(out))
  print(paste0("Read 10X : ", sam, " : ", dim(out)[1], " ", dim(out)[2]))
  test.seu <- CreateSeuratObject(counts = out, min.cells = min.cells)
  test.seu[["percent.mt"]] <- PercentageFeatureSet(test.seu, pattern = "^MT-")
  test.seu <- subset(test.seu, subset= nCount_RNA > nCount_RNA.min & percent.mt < percent.mito.max & 
                       nFeature_RNA > nFeature.min   & nCount_RNA < nCount_RNA.max &
                       nFeature_RNA < nFeature.max)
  
  test.seu <- NormalizeData(test.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  test.seu <- FindVariableFeatures(test.seu, selection.method = "vst", nfeatures = 3000)
  test.seu <- ScaleData(test.seu, features = rownames(test.seu))
  test.seu <- RunPCA(test.seu, features = VariableFeatures(test.seu),npcs = 50)
  test.seu <- FindNeighbors(test.seu, dims = 1:25)
  test.seu <- FindClusters(test.seu, resolution = 0.9)
  test.seu <- RunUMAP(test.seu, dims = 1:25)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_ovary <- paramSweep_v3(test.seu, PCs = 1:25, sct = FALSE)
  sweep.stats_ovary <- summarizeSweep(sweep.res.list_ovary, GT = FALSE)
  bcmvn_ovary <- find.pK(sweep.stats_ovary)

  ggplot(bcmvn_ovary, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()
  
  pK <- bcmvn_ovary %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- test.seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round((dim(test.seu)[2]*0.008/1000)*nrow(test.seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # run doubletFinder 
  test.seu <- doubletFinder_v3(test.seu, 
                               PCs = 1:25, 
                               pN = 0.25, 
                               pK = pK, 
                               nExp = nExp_poi.adj,
                               reuse.pANN = FALSE, sct = FALSE)
  test.seu$doublet_result <- test.seu[[colnames(test.seu@meta.data)[8]]]
  
  DimPlot(test.seu, group.by = colnames(test.seu@meta.data)[8])
  ggsave(paste0(sample.path,sam,"doublet_finder_result.png"), width = 6, height = 5, dpi = 600)
  
  Idents(test.seu) <- "doublet_result"
  VlnPlot(test.seu, c("nFeature_RNA", "nCount_RNA"))
  ggsave(paste0(sample.path,sam,"doublet_finder_vln_result.png"), width = 9, height = 5, dpi = 600)
  saveRDS(test.seu, paste0(sample.path,sam, ".rds"))
  rm(test.seu)
}

getwd()



obj.list <- stringr::str_split_fixed(sample.list, "-", 3)[,2]
ovary <- readRDS(paste0(sample.path,"OVARY-",obj.list[1],"-5P.rds"))

obj.list <- obj.list[-1]

for (sam in obj.list) {
  tmp <- readRDS(paste0(sample.path,"OVARY-",sam,"-5P.rds"))
  ovary <-merge(ovary, y = tmp)
}
dim(ovary)


sample.list <- list.files(input.path, "OV")
sample.list <- str_split_fixed(sample.list, "-", 3)[,2]

sr.raw$orig.ident <- str_split_fixed(sr.raw$orig.ident, "-", 3)[,2]
sr.raw$orig.ident <- factor(sr.raw$orig.ident, levels = rev(c("T107", "T70",
                                                                    "T04", "T110",
                                                                    "T92", "T66",
                                                                    "T17", "T18",
                                                                    "T39", "T34")))

VlnPlot(sr.raw, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident", ncol=1, pt.size = 0, combine = TRUE, cols = orig.ident.colv) & 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 20),
        axis.title.x.bottom = element_blank(),
        title = element_text(size = 23)) 
OutPutFile_QCPlot = paste0("Before_1stQC_VlnPlot_", label, ".pdf")
ggsave(file=OutPutFile_QCPlot, width=10, height=12, dpi = 600)

Idents(ovary) <- "orig.ident"

ovary$orig.ident <- factor(ovary$orig.ident, levels = rev(c("T107", "T70",
                                                              "T04", "T110",
                                                              "T92", "T66",
                                                              "T17", "T18",
                                                              "T39", "T34")))
##### after QC plot
VlnPlot(ovary, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident", ncol=1, pt.size = 0, combine = TRUE, cols = orig.ident.colv) & 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 20),
        axis.title.x.bottom = element_blank(),
        title = element_text(size = 23))
OutPutFile_QCPlot = paste0("After_1stQC_VlnPlot_", label, ".pdf")
ggsave(file=OutPutFile_QCPlot, width=10, height=12, dpi = 600)

#add primary tumor organ column
org.ids <- c("1" = "ovary", "2" = "fallopian tube", "3" = "peritoneum")

ovary$primary_tumor_organ <- ovary$orig.ident
new.organ.ids <- c(org.ids[1],org.ids[1],org.ids[1],org.ids[1],org.ids[1],org.ids[1],org.ids[2],org.ids[2],org.ids[3],org.ids[1])
ovary$primary_tumor_organ<- plyr::mapvalues(ovary$primary_tumor_organ, from = sample.list, to = new.organ.ids)

table(ovary$primary_tumor_organ)

#add specimen tissue site column

tss.ids <- c("1" = "ovary", "2" = "fallopian tube", "3" = "omentum", "4" = "lymph node, paraaortic")

ovary$specimen_tissue_site <- ovary$orig.ident
new.tissue.ids <- c(tss.ids[3],tss.ids[4],tss.ids[1],tss.ids[1],tss.ids[1],tss.ids[1],tss.ids[2],tss.ids[2],tss.ids[3],tss.ids[1])
ovary$specimen_tissue_site<- plyr::mapvalues(ovary$specimen_tissue_site, from = sample.list, to = new.tissue.ids)

table(ovary$specimen_tissue_site, ovary$orig.ident)

#cancer stage
ovary$cancer_stage <- ovary$orig.ident
new.cancer.ids <- c('3','3','4','2','2','4','1','3','4','4')
ovary$cancer_stage<- plyr::mapvalues(ovary$cancer_stage, from = sample.list, to = new.cancer.ids)

table(ovary$cancer_stage)

#metastatic
ovary$Metastatic_Primary <- ovary$orig.ident
new.meta.ids <- c('Metastatic','Metastatic','Primary','Primary','Primary','Primary','Primary','Primary','Primary','Primary')
ovary$Metastatic_Primary <- plyr::mapvalues(ovary$Metastatic_Primary, from = sample.list, to = new.meta.ids )
table(ovary$Metastatic_Primary, ovary$orig.ident)


#Recurrent/Primary

ovary$Recurrent_Primary <- ovary$orig.ident
new.meta.ids <- c('Recurrent','Recurrent','Primary','Primary','Primary','Primary','Primary','Primary','Recurrent','Primary')
ovary$Recurrent_Primary <- plyr::mapvalues(ovary$Recurrent_Primary, from = sample.list, to = new.meta.ids )

table(ovary$Recurrent_Primary)
getwd()

#Chemotherapy

ovary$chemotherapy <- ovary$orig.ident
new.meta.ids <- c('Post-treatment','Post-treatment','Post-treatment','Post-treatment',
                  'Pre-treatment','Pre-treatment','Pre-treatment','Pre-treatment',
                  'Pre-treatment','Pre-treatment')
ovary$chemotherapy <- plyr::mapvalues(ovary$chemotherapy, from = sample.list, to = new.meta.ids )


dim(ovary)
## Log Normalized data (TPM-like values)
ovary <- NormalizeData(object=ovary, normalization.method="LogNormalize", scale.factor=10000)

## Identification of highly variable features
ovary <- FindVariableFeatures(ovary, selection.method = "vst", nfeatures = num.features)
top10 <- head(VariableFeatures(ovary), 10)
dim(ovary)
#### Feature selection plot
LabelPoints(plot = VariableFeaturePlot(ovary), points = top10, repel = TRUE)
OutPutFile_QCPlot = paste0("Feature_Selection_", label, ".png")
ggsave(file=OutPutFile_QCPlot, width=6.8, height=3.8)

write.table (data.frame("VarGenes"=VariableFeatures(ovary)), file=paste0("var_genes_", label, ".txt"), row.names=FALSE, sep="\t", quote=FALSE)

## Centering the data (z-scoring)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ovary <- CellCycleScoring(ovary, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

ovary <- ScaleData(ovary, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = VariableFeatures(ovary))

## Run PCA
ovary <- RunPCA(ovary, npcs=50)

#### Elbow plot
ElbowPlot(ovary, ndims=50, reduction="pca") + theme_light()
OutPutFile_PCElbow = paste0("PCElbow_plot_dims100_", label, ".png")
ggsave(file=OutPutFile_PCElbow, width=7.5, height=6)

ovary@reductions$pca@stdev[1:41] %>%  sum() /
  ovary@reductions$pca@stdev %>% sum()
ovary@meta.data %>% View

ovary <- ovary %>% 
  RunBBKNN(run_TSNE = FALSE, n_pcs = 10, batch_key = "orig.ident") %>% 
  FindClusters(resolution = 1.5, graph.name= "bbknn") 


target.list <- list(
  Endothelial = c('PECAM1','CLDN5','FLT1','RAMP2'),
  B_Plasma = c('CD79A', 'BANK1', 'MZB1'),
  T_NKcell = c('CD3D', 'CD3E', 'TRAC','NCAM1', 'KLRD1'),
  Fibrobalst = c('DCN','BGN', 'THY1','COL1A1','COL1A2'),
  Epithelial = c('EPCAM', 'CD24', 'KRT8','KRT18', 'KRT19'),
  Myeloid = c('LYZ','CD163','CD68','FCGR3A')
)

target.features <- c('PECAM1','CLDN5','RAMP2',
                     'CD79A', 'BANK1', 'MZB1',
                     'CD3E', 'TRAC','NCAM1',
                     'DCN', 'THY1','COL1A1',
                     'EPCAM', 'KRT8','KRT18',
                     'LYZ','CD68','FCGR3A')




ovary$seurat_clusters <- factor(x = ovary$seurat_clusters, levels = c(0,17,8,23,24,
                                                                                               1,3,12,15,5,6,18,
                                                                                               2,20,
                                                                                               10,7,13,4,9,22,21,
                                                                                               11,14,19,16))
ovary$seurat_clusters <- factor(x = ovary$seurat_clusters, levels = c(0:24))

Idents(ovary) <- "seurat_clusters"
Idents(ovary) <- "seurat_clusters"


DotPlot(ovary, features = rev(paste(unlist(target.list))), cluster.idents = TRUE,
        cols= "RdYlBu",
        dot.min = 0, dot.scale = 5.5) + RotatedAxis() +
  coord_flip() +
  scale_x_discrete(position='top') +
  scale_y_discrete(position='right', limits = levels(ovary@active.ident)) +
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

ggsave("dotplot_celltype.pdf", width = 6.5, height = 6, dpi = 600)
ggsave("dotplot_celltype.png", width = 7, height = 5.5, dpi = 300)


p1 <- FeaturePlot(ovary, features= target.features, pt.size=0.01, ncol = 7) & labs(title = "") &
  theme(aspect.ratio = 1,
        legend.position='none',
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        title = element_blank())

AugmentPlot(p1)

ggsave("total_featureplot.png", width = 12, height = 6, dpi = 300)
ggsave("total_featureplot.pdf", width = 12, height = 6, dpi = 600)

##### 4. Calculated DEGs per clusters#####################################################################
##### Reload your object for appropriate PC and resolution !!!!!! ########################################

## Identify DEG

ovary$seurat_clusters <- factor(x = ovary$seurat_clusters, levels = c(0:24))
Idents(ovary) <- "seurat_clusters"

all_Markers = FindAllMarkers(object=ovary, only.pos = TRUE, min.pct = 0.25)
writexl::write_xlsx(all_Markers, path = paste0("cluster_markers_Ttest_tumor", label,".xlsx"))

all_Markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
ovary <- ScaleData(ovary, features = top10$gene)

DoHeatmap(object=ovary, features=top10$gene) + NoLegend() + theme(text = element_text(size = 5.5))
OutPutFile_Marker_heatmap = paste0("Marker_heatmap_Ttest_", label, ".png")
ggsave(file=OutPutFile_Marker_heatmap, width=33, height=15) 

##### 5. Define cell types ###############################################################################
## Write your cell type name for each cluster
celltype.id <- c("0"="Myeloid cells", "1"="Epithelial cells", "2"="Fibroblast", "3"="Epithelial cells",
                 "4"="T_NK cells", "5"="Epithelial cells", "6"="Epithelial cells", "7"="T_NK cells",
                 "8"="Myeloid cells", "9"="T_NK cells", "10"="T_NK cells","11" = "B_Plasma cells", "12"="Epithelial cells", 
                 "13"="T_NK cells", "14"="B_Plasma cells", "15"="Epithelial cells", "16"="Unknown", "17"="Myeloid cells",
                 "18"="Epithelial cells", "19" = "Endothelial cells","20" ="Fibroblast","21" ="T_NK cells","22" ="T_NK cells",
                 "23" = "Myeloid cells", "24" = "Myeloid cells")

identical(length(celltype.id), length(unique(ovary$seurat_clusters)))
ovary@meta.data$celltype <- plyr::mapvalues(ovary@meta.data$seurat_clusters, from=names(celltype.id), to=celltype.id)
ovary$celltype %>% table

celltype_table <- table(ovary$orig.ident, ovary$celltype) %>%  as.matrix.data.frame()
row.names(celltype_table) <- names(ovary$orig.ident%>% table)
colnames(celltype_table) <- names(ovary$celltype%>% table)
celltype_table <- celltype_table[,-11]

ovary$orig.ident <- str_split_fixed(ovary$orig.ident, "-", 3)[,2]



## Write your cell type order and color code
cell.type_getPalette = colorRampPalette(brewer.pal(8, "Spectral"))(10)
brewer.pal(6, "Set1")

celltype.order <- c("Epithelial cells", "Endothelial cells",
                    "Fibroblast", "T/NK cells", "B/Plasma cells", "Myeloid cells", "Unknown")

celltype.colv <- c("Epithelial cells"="#E41A1C",
                   "T/NK cells"="#377EB8",
                   "Myeloid cells"="#4DAF4A",
                   "B/Plasma cells"="#984EA3",
                   "Endothelial cells"="#FF7F00",
                   "Fibroblast"="#A65628",
                   "Unknown" = "darkgrey")

ovary$cell_lineage %>% table

orig.ident.colv <- c("T70" = "#5377AE", 
                     "T04" = "#114292", 
                     "T107" = "#5C5EAC",
                     "T110" = "#A779C5",  
                     "T66" = "#E34F79",
                     "T92" = "#FA6D3E", 
                     "T17" = "#FE9920", 
                     "T18" = "#F4CB2A", 
                     "T34" = "#ABD84C", 
                     "T39" = "#0BC7D2")

chemo_color = c("Post-treatment" = "#00B050", "Pre-treatment" = "#00B0F0")

Idents(ovary) <- "celltype"
pt.size = 0.5
library(cowplot)
library(patchwork)

#############################################################################################
p1 <- DimPlot(ovary, reduction = "umap", group.by = "cell_lineage", pt.size= pt.size, cols = celltype.colv) + 
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
ggsave("dimplot_celltype_orig_chemo.png", width = 6, height = 6, dpi = 300)
ggsave("dimplot_celltype_orig_chemo.pdf", width = 6, height = 6, dpi = 600)



legend <- get_legend(p1 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_orig_chemo_legend.png", width = 1.5, height = 3, dpi = 300)


################UMAPplot by Sample#########################
p2 <- DimPlot(ovary, reduction = "umap", group.by = "orig.ident", pt.size= pt.size, cols = orig.ident.colv) + 
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
ggsave("dimplot_celltype_orig.png", width = 6, height = 6, dpi = 300)
ggsave("dimplot_celltype_orig.pdf", width = 6, height = 6, dpi = 600)



legend <- get_legend(p2 + theme(legend.background = element_blank(),
                                legend.title = element_blank()) +  
                       guides(colour = guide_legend(override.aes = list(size=5), ncol = 1, byrow = TRUE)))
ggdraw(legend, xlim = c(1,1) )
ggsave("dimplot_celltype_orig_legend.png", width = 2, height = 6, dpi = 300)
ggsave("dimplot_celltype_orig_legend.pdf", width = 2, height = 6, dpi = 600)


p3 <- DimPlot(ovary, group.by = "seurat_clusters", reduction = "umap", pt.size= pt.size) + 
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
AugmentPlot(p3)
ggsave("total_cluster_dimplot.png", width = 6, height = 6, dpi = 300)
ggsave("total_cluster_dimplot_nolabel.pdf", width = 6, height = 6, dpi = 600)

##########UMAPplot by specimen tissue site#####################
ovary$specimen_tissue_site <- factor(ovary$specimen_tissue_site, levels = c("ovary", "omentum", "lymph node, paraaortic", "fallopian tube"))

p4 <- DimPlot(ovary, reduction = "umap", group.by = "specimen_tissue_site", pt.size=pt.size, label = FALSE, repel = TRUE, cols = c(brewer.pal(4, "Set2"))) + 
  labs(title = "") +
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

AugmentPlot(p4)
ggsave("dimplot_celltype_orig_sample_origin.png", width = 6, height = 6, dpi = 300)

legend <- get_legend(p4 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_orig_sample_origin_legend.png", width = 3, height = 1, dpi = 300)

DimPlot(ovary, group.by = "Metastatic_Primary", split.by = "orig.ident", cols = c("#e4572e", "#29335c"))
ovary$Metastatic_Primary <- factor(ovary$Metastatic_Primary, levels = c("Primary", "Metastatic"))
p4 <- DimPlot(ovary, group.by = "Metastatic_Primary", pt.size=pt.size, label = FALSE, repel = TRUE, cols = c("#e4572e", "#29335c")) + 
  labs(title = "") +
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

AugmentPlot(p4)

ggsave("dimplot_celltype_orig_sample_origin.pdf", width = 6, height = 6, dpi = 300)

legend <- get_legend(p4 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_celltype_orig_sample_origin_legend.pdf", width = 2, height = 1, dpi = 300)


p5 <- DimPlot(ovary, reduction = "umap", group.by = "doublet_result", pt.size=pt.size, label = FALSE, repel = TRUE, cols = c("lightgrey","red")) + 
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

AugmentPlot(p5)
ggsave("dimplot_doublet.png", width = 6, height = 6, dpi = 300)

legend <- get_legend(p5 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_doublet_legend.png", width = 2, height = 1, dpi = 300)

ovary$chemotherapy %>% table
ovary$chemotherapy <- factor(ovary$chemotherapy, levels = c("Pre-treatment", "Post-treatment"))
p6 <- DimPlot(ovary, reduction = "umap", group.by = "chemotherapy", pt.size=pt.size, cols = c("#00B0F0", "#00B050")) + 
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
ggsave("dimplot_total_treatment.pdf", width = 6, height = 6, dpi = 600)

legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE)))
ggdraw(legend)
ggsave("dimplot_total_tratment_legend.pdf", width = 4, height = 1, dpi = 600)

#####################################################################
######### Percent bar graph (by orig.ident - clusters) ##############
library(dplyr)
library(RColorBrewer)
library(gghighlight)
library(tidyr)
library(ggsignif)
library(gdata)

## Percent bar graph (by cell.type - orig.ident)
Input.df <- ovary@meta.data
nCells2 <- dplyr::count (Input.df, orig.ident)
pre_Counted2 <- Input.df %>% group_by(orig.ident) %>% dplyr::count(cell_lineage)
Counted2 <- pre_Counted2 %>% ungroup %>% tidyr::complete(cell_lineage, orig.ident,  fill = list(n = 0)) # with zero-values

Input2 <- merge (x=Counted2, y=nCells2, by = "orig.ident", all.x=TRUE)
Input2$Percent <- Input2$n.x / Input2$n.y * 100
table(Input2$cell_lineage)


Input2$orig.ident <- factor(Input2$orig.ident, levels = c("T107", "T70",
                                                          "T04", "T110",
                                                          "T92", "T66",
                                                          "T17", "T18",
                                                          "T39", "T34"))

count_input <- unique(Input2[,c(1,4)])
count_input$orig.ident <- factor(count_input$orig.ident, levels = c("T107", "T70",
                                                                    "T04", "T110",
                                                                    "T92", "T66",
                                                                    "T17", "T18",
                                                                    "T39", "T34"))



Input2$cell_lineage <- factor(Input2$cell_lineage, levels=(rev(c("Unknown", "Fibroblast", "Endothelial cells", 
                                                         "B/Plasma cells","Myeloid cells", "T/NK cells", "Epithelial cells"))))


p1 <- ggplot(Input2, aes(x=orig.ident, y=Percent, fill=cell_lineage)) + 
  geom_bar(stat="identity",position = position_fill(reverse = TRUE), color="black")  +
  labs(x=NULL, y="% Proportion") + coord_flip() + 
  theme_classic() + theme(axis.text.y = element_text(size = 15),
                          axis.text.x = element_text(size = 11), 
                          axis.title.x = element_text(size = 13),
                          legend.position = "bottom",
                          legend.text = element_text(size =13)) +
  scale_fill_manual(values = celltype.colv) +   
  guides(fill = guide_legend(override.aes = list(size=3), ncol = 2, byrow = TRUE))

p2 <- ggplot(count_input, aes(x=orig.ident, y=n.y)) + geom_bar(stat="identity")  + 
  labs(x=NULL, y="Count") + coord_flip() + NoLegend() +
  theme_classic() + theme(aspect.ratio= 2.3,
                          axis.ticks = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size = 11), 
                          axis.title.x = element_text(size = 13))

p1 + p2

OutPutFile_bar = paste("total_Clusters_Samples_orig.ident_legend.png", sep="")
ggsave(file=OutPutFile_bar, width =7, height = 6 ,dpi=300)
OutPutFile_bar = paste("total_Clusters_Samples_orig.ident_legend.pdf", sep="")
ggsave(file=OutPutFile_bar, width =7, height = 6 ,dpi=600)


getwd()

saveRDS(ovary, "ovary_total.RDS")

library(gdata)
Idents(ovary) <- "orig.ident"

pair <- subset(ovary, idents = c("T92", "T110"))
table(pair$celltype)
Idents(pair) <- "celltype"
pair <- subset(pair, idents = c("Myeloid cells", "Fibroblast", "T_NK cells", "B_Plasma cells", "Endothelial cells"))

pair$celltype <- drop.levels(pair$celltype)
table(pair$celltype) %>% names()


for (celltype in table(pair$celltype) %>% names()) {
  lineage <- subset(pair, idents = celltype)
  Idents(lineage) <- "orig.ident"
  clustering <- lineage@active.ident
  cellsIn <- names(clustering[clustering == "T110"])
  cellsOut <- names(clustering[clustering == "T92"])
  lineage<- ScaleData(lineage)
  
  markers <- FindMarkers(lineage, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)
  markers$genes <- rownames(markers)
  writexl::write_xlsx(markers, path = paste0("project_directory_name/ovary_project/Seurat/final/", "DEG",celltype, "pair_sample",".xlsx"))
  endo_markers$p_val_adj
  EnhancedVolcano(markers, lab = rownames(markers), x = 'avg_log2FC', y = 'p_val_adj', title = paste0("T110 vs T92" , celltype))
  ggsave(paste0("project_directory_name/ovary_project/Seurat/final/","T110 vs T92" , celltype,".png"), width = 15, height = 12)
}

table(pair$orig.ident, pair$celltype)

lineage <- subset(pair, idents = "Fibroblast")
Idents(lineage) <- "orig.ident"
clustering <- lineage@active.ident
cellsIn <- names(clustering[clustering == "T110"])
cellsOut <- names(clustering[clustering == "T92"])
lineage<- ScaleData(lineage)

markers <- FindMarkers(lineage, ident.1 = cellsIn, ident.2 = cellsOut, test.use = 'MAST', logfc.threshold = 0)


lineage$orig.ident %>% table
DimPlot(lineage)
p1 <- DimPlot(fibro, group.by = "celltype") 
p2 <- FeaturePlot(fibro, "C7", split.by = "orig.ident") 
p2

pair_down_genes<- c("ECM1", "FXYD5", "GPNMB", "CAPG", "EMP1", "EGFL6", "CRABP2", "AEBP1", "GAPDH", "CD99", "TUBA1C", "IGF1", "COL8A1", "GJB2", "OLFML2B", "OSTC", "UPK3B", "MMP14", "THY1", "VCAN", "S100A6", "S100A11", "SH3BGRL3", "FBLN1", "CXCL14", "DES", "ANXA1", "DERL3", "S100A10", "TUBA1B", "HP", "PLIN2", "C3", "THBS2", "COL5A1", "PTGDS", "PLAT", "CTSK", "CALB2", "KRT19", "H19", "ANXA2", "TUBA1A", "ACTG2", "COL1A2", "FN1", "FAP", "SULF1", "CTHRC1", "LOX", "SPARC", "MFAP5", "C1QTNF3", "COL5A2", "KRT8", "CFB", "PLAU", "MFAP2", "COL1A1", "COL3A1", "SERPINB2", "COL11A1", "POSTN", "LUM", "EPYC", "CCDC80", "SLPI", "SFRP2", "MMP11")

total_down_genes <- c("MMP11", "ISG15", "TIMP1", "IFIT1", "IFI6", "EGFL6", "IFI27", "COL4A1", "HLA-A", "HLA-C", "SLPI", "B2M", "COL4A2", "TDO2", "CXCL10", "EPYC", "PTGDS", "ACTG2", "IFIT3", "IFI44L", "PRSS23", "BST2", "COL18A1", "STAT1", "DES", "SH3BGRL3", "MX1", "IFITM1", "CYTOR", "NDUFA4L2", "MYL12A", "DERL3", "PSME2", "STMN1", "LY6E", "CFB", "GAPDH", "OAS1", "TMSB4X", "XAF1", "MMP19", "WFDC2", "IFI30", "PFN1", "BNIP3", "CCND1", "MIF", "THY1", "PLIN2", "POMP", "PLAT", "LGALS3BP", "HSPA6", "IFI44", "DUXAP8", "APOC1", "SSR3", "S100A11", "DYNLT1", "EPSTI1", "CKS2", "ITGA1", "CD74", "ARPC3", "ADM", "IRF7", "KRT19", "HSP90B1", "VAMP5", "H2AFZ", "CHN1", "IGFBP2", "MX2", "CALD1", "OAS3", "PARP14", "TNFAIP6", "COLEC11", "VMP1", "COL6A1", "MYLK", "MALAT1", "XIST", "RPS27L", "RNF213", "TMSB10", "HSPB1", "CALM2", "GBP1", "PPIB", "OASL", "VEGFA", "CFL1", "CA12", "GPM6B", "PGK1", "CHCHD10", "IFIT2", "PLSCR1", "ECM1", "CAV1", "SERPINH1", "TPI1", "COL5A3", "RNASE1", "CLIC1", "NEAT1")

pair_down_genes[pair_down_genes%in% total_down_genes]
total_down_genes[total_down_genes %in%pair_down_genes]

pair_up_genes<-c("STAR", "COLEC11", "MT-ATP8", "C7", "RBP1", "RGS5", "MT-ATP6", "DDIT4", "TSC22D1", "MT-ND4L", "MT-CO1", "TSPAN8", "IGFBP3", "CEBPD", "MT-ND1", "DNAJB1", "CCL2", "DEPP1", "FHL2", "CXCL8", "FXYD6", "PPP1R14A", "NR4A1", "MT-CO3", "MTRNR2L12", "CCN1", "GJA4", "ARHGAP29", "RNASE1", "MT-CYB", "MT-ND5", "DKK3", "JUND", "PEG3", "PDLIM3", "LITAF", "RHOB", "AKAP12", "ID2", "ADIRF", "SLC40A1", "IGFBP7", "LGR5", "SOCS3", "MT-ND3", "ZFP36L2", "SGK1", "RGS16", "MT-ND4", "BCAM", "TCF21", "SERPINA3", "GADD45B", "MCAM", "ADAMTS4", "SERPINA5", "CXCL3")

total_up_genes <-c("SFRP2", "CFD", "MGP", "C7", "MFAP4", "CXCL14", "SFRP4", "ADH1B", "FBLN1", "PDK4", "GSN", "MTRNR2L12", "IGFBP5", "C3", "ADIRF", "ELN", "APOD", "DPT", "TXNIP", "KLF4", "CCN2", "BGN", "IGFBP4", "PODN", "MT-ATP6", "MT1M", "IGFL2", "TSC22D1", "TIMP3", "CCN5", "IER2", "MT-ND3", "CYP1B1", "SEMA3C", "MT-ND4", "HSPB6", "CLU", "MT-ND4L", "MT-CYB", "MT-ND5", "ALDH2", "FBLN2", "COMP", "CCN1", "JUN", "ASPN", "EGR1", "COL14A1", "GADD45B", "GAS6", "TSC22D3", "MT1A", "MT1X", "DNAJB1", "PCSK1N", "MT-ND1", "CCDC80", "LRP1", "AHNAK", "BOC", "MMP2", "SERPINF1", "MT-CO3", "RPS3A", "CRISPLD2", "FOS", "CLDN11", "COL8A1", "COL16A1", "PLAC9", "DEPP1", "IGFBP3", "FBLN5", "ZBTB16", "FOSB", "CHRDL1", "FIBIN", "CILP", "ZFP36L2", "MT-CO1", "JUNB", "PLK2", "PTGER3", "MYC", "CEBPD", "LAMP5", "FXYD1", "SSPN", "SPON1", "CCL2", "MT-ND2", "SOCS3", "COL3A1", "PDGFRL", "FMO2", "TMEM119", "MT2A", "PPP1R15A", "RARRES1", "DCN", "SRPX", "RPS27", "CD34", "NOP53", "PRELP", "RPL23A", "RPL37A", "IL6", "NFIB", "LTBP4", "PLA2G2A", "LOXL1", "GLUL", "GPNMB", "ZFP36", "RPS2")

pair_up_genes[pair_up_genes%in% total_up_genes]
total_up_genes[total_up_genes %in%pair_up_genes]

FeaturePlot(endo, "CD36", split.by = "orig.ident")

