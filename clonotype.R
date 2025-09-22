# R-4.0.2
remove(list = ls())

library(scRepertoire) # v1.1.2 # devtools::install_github("ncborcherding/scRepertoire")
library(Seurat)
library(RColorBrewer)
library(ggalluvial)
library(ggplot2)
library(dplyr)

# Control Part //////////////////////////////////////////////////////#
dir <- "D:/project directory name/scanpy/ovarian/tnk/TCR/"
ipath <- paste0(dir)
opath <- paste0(dir)


tnk_obj_path <- paste0(dir, "../")
tnk_obj <- "1.3celltype_annotation_Tcell.rds"

#////////////////////////////////////////////////////////////////////#

setwd(opath)

# 'Frequency' would be calculated based on cloneCall/CTtype
CTtype <- "CTstrict"
cloneCall <- "gene+nt"

# load seurat obj
tnk <- readRDS(paste0(tnk_obj_path, tnk_obj))
tnk <- SetIdent(tnk, value = "celltype")
table(Idents(tnk))
tnk@meta.data$celltype %>%  table
tnk <- readRDS("d:/project directory name/scanpy/ovarian/tnk/TCR/CombineExpression_CTstrict.rds")

# load scReperotire obj
combined_filter <- readRDS("scRepertoire/CombineTCR.rds")
combined_filter2 <- readRDS("CombineTCR_OnlyPairs.rds")

View(combined_filter[[1]])
### 6. Interacting with Seurat #######################################
# renaming labels
combined_filter[1]
for (i in 1: length(combined_filter)){
  combined_filter[[i]]$barcode <- gsub("__", "_", combined_filter[[i]]$barcode)
  
}
View(combined_filter[[1]])

# Combine to expression data
# -> 'Frequency' calculation based on cloneCall type !!!
seurat <- combineExpression(combined_filter, tnk, cloneCall=cloneCall, groupBy = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=10, Large=20, Hyperexpanded=Inf))
tnk$orig.ident %>% table
table(seurat$cloneType)

# Assign cloneType by frequency ranges
slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (20 < X <= Inf)", "Large (10 < X <= 20)", 
                                                         "Medium (5 < X <= 10)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))
View(seurat@meta.data)

# save RDS ----------------------------------------------------------
saveRDS(seurat, file = paste0("CombineExpression_", CTtype, ".rds"))
#--------------------------------------------------------------------

## UMAP by cloneType ---##
# color key
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
clonetype.col <- colorblind_vector(5)
names(clonetype.col) <- levels(seurat$cloneType)
clonetype.col["NA"] <- "gray"

DimPlot(seurat, group.by = "cloneType") + scale_color_manual(values = clonetype.col, na.value = "lightgray") +#scale_color_manual(values = colorblind_vector(5), na.value="grey")
  theme(axis.title = element_text(size = 10),
        axis.text=element_blank(), 
        axis.ticks=element_blank()) + labs(x = "bbknnumap_1", y = "bbknnumnap_2")
  ggsave(paste0("Seurat_UMAP_cloneType_", CTtype, "_legend.png"), dpi=300, width = 8, height = 4)
DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = clonetype.col, na.value = "lightgray") + #scale_color_manual(values = colorblind_vector(5), na.value="grey")
  NoLegend()
ggsave(paste0("Seurat_UMAP_cloneType_", CTtype, ".png"), dpi=300, width = 4, height = 4)


##samples barplot##
tnk <- SetIdent(tnk, value = 'highlight_top5')
meta <- tnk@meta.data
meta <- tnk[!is.na(meta$cloneType),]

ncells <- count(meta, orig.ident)
pre_count <- meta %>% group_by(orig.ident) %>% count(cloneType)
count <- merge(pre_count,ncells, by = "orig.ident")
count$prop <- count$n.x/count$n.y *100


count$orig.ident <- factor(count$orig.ident, levels = c("OVARY-T70-5P", "OVARY-T04-5P", 
                                                        "OVARY-T107-5P", "OVARY-T110-5P",  
                                                        "OVARY-T92-5P", "OVARY-T66-5P",  
                                                        "OVARY-T17-5P", "OVARY-T18-5P", 
                                                        "OVARY-T34-5P", "OVARY-T39-5P"))
count2 <- unique(count[,c(1,4)])


p3 <- ggplot(count, aes(fill = cloneType, x = orig.ident, y = prop )) + coord_flip() +
  geom_bar(color="black",stat="identity", position = position_fill(reverse = TRUE)) + 
  labs(x=NULL, y="% proportion") +
  theme_classic() + theme(legend.position = "bottom", 
                          legend.key.width = unit(0.4,"cm"),
                          legend.key.size =  unit(0.4,"cm"),
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(colour="black", size= 7)) +
  scale_fill_manual(values = c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")) 

p3 + ggplot(count2, aes(x = orig.ident, y = n.y )) + coord_flip() +
  geom_bar(color="black",stat="identity") + labs(x=NULL, y="count") +
  theme_classic() + theme(aspect.ratio= 3,
                          axis.ticks = element_blank(),
                          legend.position = "right", 
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          axis.text.y = element_blank())

OutPutFile_bar = paste("t_clonotype_sample.png", sep="")
ggsave(file=OutPutFile_bar, width = 8, height = 6   ,dpi=600)

######celltype#####

meta <- seurat@meta.data
meta <- meta[!is.na(meta$cloneType),]

ncells <- count(meta, celltype)
pre_count <- meta %>% group_by(celltype) %>% count(cloneType)
count <- merge(pre_count,ncells, by = "celltype")
count$prop <- count$n.x/count$n.y *100
count2 <- count[,c(1,4)] %>% unique()

p4 <- ggplot(count, aes(fill = cloneType, x = celltype, y = prop )) + coord_flip() +
  geom_bar(color="black",stat="identity", position = position_fill(reverse = TRUE)) + 
  labs(x=NULL, y="%proportion") +
  theme_classic() + theme(legend.position = "bottom", 
                          legend.key.width = unit(0.4,"cm"),
                          legend.key.size =  unit(0.4,"cm"),
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(colour="black", size= 7)) +
  scale_fill_manual(values = c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")) +
  scale_x_discrete(limit=c("CD8+ GZMB",
                           "CD8+ GZMK",
                           "CD8+ FGFBP2",
                           "CD8+ IL7R",
                           "CD4+ CXCL13",
                           "CD4+ FOXP3",
                           "CD4+ IL7R",
                           "NK cells",
                           "Proliferating"))


p4 + ggplot(count2, aes(x = celltype, y = n.y )) + coord_flip() +
  geom_bar(color="black",stat="identity") + labs(x=NULL, y="count") +
  theme_classic() + theme(aspect.ratio= 3,
                          axis.ticks = element_blank(),
                          legend.position = "right", 
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          axis.text.y = element_blank()) +
  scale_x_discrete(limit=c("CD8+ GZMB",
                           "CD8+ GZMK",
                           "CD8+ FGFBP2",
                           "CD8+ IL7R",
                           "CD4+ CXCL13",
                           "CD4+ FOXP3",
                           "CD4+ IL7R",
                           "NK cells",
                           "Proliferating"))

OutPutFile_bar = paste("t_clonotype_celltype.png", sep="")
ggsave(file=OutPutFile_bar, width = 8, height = 6   ,dpi=600)

####boxplot####


