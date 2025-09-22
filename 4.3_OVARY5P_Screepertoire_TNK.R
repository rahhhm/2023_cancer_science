# R-4.0.2

library(scRepertoire) # v1.1.2 # devtools::install_github("ncborcherding/scRepertoire")
library(Seurat)
library(RColorBrewer)
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(cowplot)

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
combined_filter <- readRDS("CombineTCR.rds")
combined_filter2 <- readRDS("CombineTCR_OnlyPairs.rds")

View(combined_filter[[1]])
### 6. Interacting with Seurat #######################################

# Combine to expression data
# -> 'Frequency' calculation based on cloneCall type !!!

tnk <- combineExpression(combined, tnk, cloneCall=cloneCall, group.by = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=10, Large=Inf))

# Assign cloneType by frequency ranges
slot(tnk, "meta.data")$cloneType <- factor(slot(tnk, "meta.data")$cloneType, 
                                              levels = c("Large (10 < X <= Inf)", 
                                                         "Medium (5 < X <= 10)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))

tnk$cloneType %>% table()

# save RDS ----------------------------------------------------------
saveRDS(seurat, file = paste0("CombineExpression_", CTtype, ".rds"))
#--------------------------------------------------------------------

## UMAP by cloneType ---##
# color key
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF"))
clonetype.col <- colorblind_vector(4)
names(clonetype.col) <- levels(tnk$cloneType)
clonetype.col["NA"] <- "gray"

tnk$cloneType <- factor(tnk$cloneType, levels = c("Large (10 < X <= Inf)", 
                                         "Medium (5 < X <= 10)", "Small (1 < X <= 5)", "Single (0 < X <= 1)"))

tnk$chemo  <- factor(tnk$chemo, levels = c("Pre-treatment", "Post-treatment"))
p6 <- DimPlot(tnk, group.by = "cloneType", pt.size = 1, split.by = "chemo")  +
  labs(title = "", subtitle = "", caption = "", tag = "", alt = "", alt_insight = "") +
  scale_color_manual(values = clonetype.col, na.value = "lightgray") + 
  theme(aspect.ratio= 1,
        line = element_blank(), 
        rect = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_blank()) 

AugmentPlot(p6)

ggsave(paste0("clonotype_splited", CTtype, ".png"), width = 20, height = 20, dpi = 300)
ggsave(paste0("clonotype_splited", CTtype, ".pdf"), width = 11, height = 3.5, dpi = 600)


legend <- get_legend(p6 + theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.text = element_text(colour="black", size=10)) +  
                       guides(colour = guide_legend(override.aes = list(size=3), ncol = 4, byrow = TRUE)))
ggdraw(legend)
ggsave("tnk+dimplot_total_treatment_legend.pdf", width = 8, height = 1, dpi = 600)


# save RDS ----------------------------------------------------------
saveRDS(seurat, file = paste0("CombineExpression_", CTtype, ".rds"))
#--------------------------------------------------------------------

##samples barplot##
  
grep("NA",meta$CTnt)
meta <- tnk@meta.data 
meta <- meta[!is.na(meta$cloneType),]
meta <- meta %>% filter(celltype != "NK" & orig.ident != "T39")
meta$celltype %>% table
ncells <- dplyr::count(meta, orig.ident)
pre_count <- meta %>% group_by(orig.ident) %>% dplyr::count(cloneType)
count <- merge(pre_count,ncells, by = "orig.ident")
count$prop <- count$n.x/count$n.y *100


count$orig.ident <- factor(count$orig.ident, levels = c("T107", "T70",
                                                        "T04", "T110",
                                                        "T92", "T66",
                                                        "T17", "T18",
                                                        "T34"))
count2 <- unique(count[,c(1,4)])

p3 <- ggplot(count, aes(fill = cloneType, x = orig.ident, y = prop )) + coord_flip() +
  geom_bar(color="black",stat="identity", position = position_fill(reverse = TRUE)) + 
  labs(x=NULL, y="% Proportion") +
  theme_classic() + theme(legend.position = "bottom",
                          legend.direction = "horizontal",
                          legend.key.width = unit(0.4,"cm"),
                          legend.key.size =  unit(0.4,"cm"),
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          legend.text = element_text(colour="black", size= 19),
                          axis.text = element_text(size = 15)) +
  guides(fill=guide_legend(ncol=2,byrow=TRUE)) +
  scale_fill_manual(values = c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF")) 

p3 + ggplot(count2, aes(x = orig.ident, y = n.y )) + coord_flip() +
  geom_bar(color="black",stat="identity") + labs(x=NULL, y="Count") +
  theme_classic() + theme(aspect.ratio= 3,
                          axis.ticks = element_blank(),
                          legend.position = "right", 
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          axis.text.y = element_blank(),,
                          axis.text = element_text(size = 15))

OutPutFile_bar = paste("t_clonotype_sample.pdf", sep="")
ggsave(file=OutPutFile_bar, width = 8, height = 6   ,dpi=600)

######celltype#####

meta <- tnk@meta.data
meta <- meta %>% filter(celltype != "NK" & orig.ident != "T39")
meta$celltype <- gdata::drop.levels(meta$celltype)
meta$orig.ident <- gdata::drop.levels(meta$orig.ident)

meta$celltype <- factor(meta$celltype, levels = c("CD4 Tn", "CD4 T ISG", "CD8 Tn/m", "CD4 Tm", "CD4 Treg", "CD8 T ISG", "CD8 Tem", "CD4 Tfh","Proliferative cells","CD8 Tex"))


meta <- meta[!is.na(meta$cloneType),]
meta$celltype %>% table

ncells <- dplyr::count(meta, celltype)
pre_count <- meta %>% group_by(celltype) %>% dplyr::count(cloneType)
count <- merge(pre_count,ncells, by = "celltype")
count <- na.omit(count)
count$prop <- count$n.x/count$n.y *100
count2 <- count[,c(1,4)] %>% unique()

p4 <- ggplot(count, aes(fill = cloneType, x = celltype, y = prop )) + coord_flip() +
  geom_bar(color="black",stat="identity", position = position_fill(reverse = TRUE)) +
  theme_classic() + labs(x=NULL, y="% Proportion") + theme(legend.position = "bottom",
                                         legend.direction = "horizontal",
                                         legend.key.width = unit(0.4,"cm"),
                                         legend.key.size =  unit(0.4,"cm"),
                                         legend.background = element_blank(),
                                         legend.title = element_blank(),
                                         legend.text = element_text(colour="black", size= 16),
                                         axis.text = element_text(size = 12)) +
  guides(fill=guide_legend(ncol=2,byrow=TRUE)) +
  scale_fill_manual(values = c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF"))

p4 + ggplot(count2, aes(x = celltype, y = n.y )) + coord_flip() +
  geom_bar(color="black",stat="identity") + labs(x=NULL, y="Count") +
  theme_classic() + theme(aspect.ratio= 3,
                          axis.ticks = element_blank(),
                          legend.position = "right", 
                          legend.background = element_blank(),
                          legend.title = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text = element_text(size = 12))

OutPutFile_bar = paste("t_clonotype_celltype.pdf", sep="")
ggsave(file=OutPutFile_bar, width = 9, height = 6, dpi=300)
ggsave(file=OutPutFile_bar, width = 9, height = 6, dpi=600)


combined2 <- expression2List(tnk_pair, group = "celltype")

length(combined2) #now listed by cluster
Idents(tnk) <- "orig.ident"
tnk_pair <- subset(tnk, idents = c("T92", "T110"))

clonalOverlap(combined2, 
              cloneCall="aa", 
              method="overlap") + rotate_x_text()
ggsave("overlapped_clonotype_seq_celltype_pair.png", width = 10, height = 8, dpi = 300)

clonalOverlap(combined2, 
              cloneCall = "aa", 
              method = "morisita", chain = "both") 

ggsave("overlapped_clonotype_seq_celltype_pair.png", width = 6.5, height = 5, dpi = 300)


####Clonotype boxplot####
library(tidyverse)
library(ggpubr)

meta <- tnk@meta.data 
meta <- meta[!is.na(meta$cloneType),]
meta$clonotype <- meta$cloneType
ncells <- dplyr::count(meta, orig.ident)
pre_count <- meta %>% group_by(orig.ident) %>% dplyr::count(clonotype, .drop = FALSE)
count <- merge(pre_count,ncells, by = "orig.ident")
count$prop <- count$n.x/count$n.y *100


count$orig.ident <- factor(count$orig.ident, levels = rev(c("T107", "T70",
                                                        "T04", "T110",
                                                        "T92", "T66",
                                                        "T17", "T18",
                                                        "T39", "T34")))

count$orig.ident

sample.list <- rev(c("T107", "T70",
                 "T04", "T110",
                 "T92", "T66",
                 "T17", "T18",
                 "T39", "T34"))

count$chemo <- count$orig.ident
new.meta.ids <- c('Pre-treatment','Pre-treatment','Pre-treatment','Pre-treatment',
                  'Pre-treatment','Pre-treatment', 'Post-treatment','Post-treatment',
                  'Post-treatment', 'Post-treatment')
count$chemo <- plyr::mapvalues(count$chemo, from = sample.list, to = new.meta.ids )

count <- count %>% filter(orig.ident != "T39")

chemo_color <- c("#00B0F0", "#00B050")
my_comparisons = list(c("Pre-treatment", "Post-treatment"))

p1 <- count %>% filter(clonotype == "Small (1 < X <= 5)") %>%
  ggplot(aes(x = chemo, prop, fill = chemo)) + theme_bw() +
  geom_boxplot(width = 0.5) +
  geom_point(size = 1.5, shape = 19) + geom_jitter(width = 0, color = "black") +
  labs(title = "Small", y = "")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.title.y.left = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        title = element_text(size = 13)) +
  scale_fill_manual(values = chemo_color) +  NoLegend() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4)

p2 <- count %>% filter(clonotype == "Medium (5 < X <= 10)") %>%
  ggplot(aes(x = chemo, prop, fill = chemo)) + theme_bw() +
  geom_boxplot(width = 0.5) +
  geom_point(size = 1.5, shape = 19) + geom_jitter(width = 0, color = "black") +
  labs(title = "Medium", y = "")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.title.y.left = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        title = element_text(size = 13)) +
  scale_fill_manual(values = chemo_color) + NoLegend() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4)

p3 <- count %>% filter(clonotype == "Large (10 < X <= Inf)") %>%
  ggplot(aes(x = chemo, prop, fill = chemo)) + theme_bw() +
  geom_boxplot(width = 0.5) +
  geom_point(size = 1.5, shape = 19) + geom_jitter(width = 0, color = "black") +
  labs(title = "Large", y = "")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.title.y.left = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        title = element_text(size = 13)) +
  scale_fill_manual(values = chemo_color) + NoLegend() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4)

p4 <- count %>% filter(clonotype == "Single (0 < X <= 1)") %>%
  ggplot(aes(x = chemo, prop, fill = chemo)) + theme_bw() +
  geom_boxplot(width = 0.5) +
  geom_point(size = 1.5, shape = 19) + geom_jitter(width = 0, color = "black") +
  labs(title = "Single", y = "% Proportions")+
  theme(panel.background = element_rect(inherit.blank = F),
        axis.title.y.left = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_blank(),
        title = element_text(size = 13)) +
  scale_fill_manual(values = chemo_color) + NoLegend() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", hide.ns = TRUE, size = 4)

p4 + p1 + p2 + p3 + patchwork::plot_layout(ncol = 4)
ggsave("clonotype_boxplot.png", width = 8, height = 4.3, dpi = 300)
ggsave("../clonotype_boxplot.pdf", width = 8, height = 5, dpi = 600)

################ 7. Overlapping rate of clonotypes between Types ###################################
# manual analysis after 'compareClonotypes' analysis
library(ggvenn)
library(dplyr)

########################### Generate filtered data##########################
meta <- tnk@meta.data
meta <- meta %>% filter(celltype != "NK cells")
f1 <- meta %>% filter(orig.ident == "T92") %>% filter(Frequency > 0)
f2 <- meta %>% filter(orig.ident == "T110") %>% filter(Frequency > 0)

########################## Generate data sets###################################
A <- list('T92 (Pre)' = c(f1$CTaa %>% unique()), 'T110 (Post)' = c(f2$CTaa %>% unique()))

ggvenn(A,show_percentage = FALSE, fill_color = c("#A779C5", "#FA6D3E"), stroke_color = "black", stroke_size = 0.1, text_size = 8, set_name_size = 8)

ggsave("venndiergram.png", width = 8, height = 7, dpi = 300)
ggsave("venndiergram.pdf", width = 8, height = 7, dpi = 600)



meta <- meta %>% filter(orig.ident == "T92" |orig.ident == "T110")

meta92 <- meta %>% filter(orig.ident == "T92")
meta110 <- meta %>% filter(orig.ident == "T110")
meta92$overlapped <- ifelse(meta92$CTaa %in% meta110$CTaa, "Overlapped", "Unique")
meta110$overlapped <- ifelse(meta110$CTaa %in% meta92$CTaa, "Overlapped", "Unique")

meta_overlapped <- rbind(meta110, meta92)

meta_overlapped <- meta_overlapped %>% filter(overlapped == "Overlapped")


meta_overlapped %>% colnames
meta_overlapped$celltype2 <- meta_overlapped$celltype

celltype.id2 = c("CD8 Tex + ISG", "NK", "CD4 Tm", "CD8 Tem",
                 "CD8 Tex + ISG", "CD4 Treg", "CD4 Tn",  "CD4 Tfh", "Proliferating", 
                 "CD8 Tn/m","CD4 T ISG")
meta_overlapped$celltype2  <- plyr::mapvalues(meta_overlapped$celltype2 , from=names(table(meta_overlapped$celltype)), to=celltype.id2)
meta_overlapped$celltype2 %>% table

meta_overlapped <- meta_overlapped %>% mutate(celltype_orig.ident = paste0(celltype2 , "_", orig.ident))
meta_overlapped %>% colnames
meta_overlapped <- table(meta_overlapped$CTaa, meta_overlapped$celltype_orig.ident) %>% as.data.frame()
meta_overlapped <- meta_overlapped %>% filter(Freq != 0)

meta_overlapped$celltype <- stringr::str_split_fixed(meta_overlapped$Var2, pattern = "_", n = 2)[,1]
meta_overlapped$orig.ident <- stringr::str_split_fixed(meta_overlapped$Var2, pattern = "_", n = 2)[,2]
meta_overlapped <- meta_overlapped[,c(1,5,4,3)]
meta_overlapped %>% colnames

meta_overlapped$celltype %>% table
meta_overlapped <- meta_overlapped %>%  filter(celltype != "NK")
meta_overlapped$celltype <- factor(meta_overlapped$celltype, levels = rev(c("Proliferating", "CD4 Tm",  "CD8 Tn/m",  "CD4 T ISG", "CD8 Tem",
                                                                            "CD8 Tex + ISG", "CD4 Treg")))

clonotype_list <- meta_overlapped$Var1  %>% unique
clonotype_list <- clonotype_list[rev(c(6,7,8,9,21,22,19,15,16,10,12,14,18,23,17,13,11,5,20,1,2,4,3))]
# p1 <- meta_overlapped %>% filter(orig.ident == "T92") %>% ggplot(aes(Var1, Freq, fill=celltype))+geom_bar(stat='identity', color = "white") + theme_classic() + rotate_x_text() + labs(x = "", y = "T92 (Pre)") + 
#   scale_fill_manual(values = celltype.colv) +  coord_flip() + scale_y_reverse(limits=c(20, 0)) + 
#   theme(aspect.ratio = 1, legend.position = "top", axis.text.y = element_text(size = 7)) + NoLegend() +
#   geom_hline(yintercept=c(0:20), linetype='solid', color='white', size=0.5) +
#   xlim(clonotype_list)
# 
# 
# p2 <-  meta_overlapped %>% filter(orig.ident == "T110") %>% ggplot( aes(Var1, Freq, fill=celltype))+geom_bar(stat='identity', color = "white") + theme_classic() + rotate_x_text() + labs(title = " ", y = "T110 (Post)") + 
#   scale_fill_manual(values = celltype.colv) +  coord_flip() +  theme(aspect.ratio = 1, axis.line.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+ NoLegend() + ylim(c(0,20)) +
#   geom_hline(yintercept=c(0:20), linetype='solid', color='white', size=0.5)+
#   xlim(clonotype_list)
# 
# l <- cowplot::get_legend(p1 + theme(legend.position = "bottom",
#                                     legend.key.size = unit(0.1, 'cm'), 
#                                     legend.key.height = unit(0.1, 'cm'),
#                                     legend.key.width = unit(0.3, 'cm')))
# 
# cowplot::plot_grid(l) / (p1 + p2 ) + 
#   plot_layout(heights = unit(c(1, 5), "cm"))
# 
# ggsave("overlapped_clonotype_barplot.png", width = 8, height = 4, dpi = 300)


celltype.colv2 <- c("CD4 Tm" = "#00b10f",
                   "CD4 T ISG" = "#e466b3",
                   "CD4 Tfh" = "#bc6c25",
                   "CD4 Treg" = "#6d2727",
                   "CD8 Tem" = "#0BC7D2",
                   "CD8 Tex + ISG" = "#e92244",
                   "Proliferating" = "#2a9d8f")


ggplot(meta_overlapped, aes(x = orig.ident, y = Var1)) + geom_point(aes(size = Freq, color = celltype) ) + 
  facet_grid(~celltype) + labs(x = "", y = "Clonotype") + theme_bw() + scale_size(range = c(2,10)) +
  scale_color_manual(values = celltype.colv2) + scale_x_discrete(limits=c("T92", "T110"))  +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ylim(clonotype_list)

ggsave("d:/project_directory_name/ovary_project/Seurat/final/tnk/facet_celltype_clonotype_dotplot.pdf", width = 15, height = 6, dpi = 600)



prop <- table(meta_overlapped$orig.ident, meta_overlapped$celltype) %>% as.data.frame()
colnames(prop) <- c("orig.ident", "celltype", "Freq")
nCells<- aggregate(prop$Freq, list(prop$orig.ident), FUN=sum) %>% as.data.frame()
colnames(nCells) <- c("orig.ident", "n")
prop <- merge(prop, nCells, by = "orig.ident")
prop$Percent <- prop$Freq / prop$n * 100
prop$orig.ident <- factor(prop$orig.ident, levels = c("T92", "T110"))
ggplot(prop, aes(x="", y=Percent, fill=celltype))+ theme_bw() + labs(title = "", x = "") +
  geom_bar(width = 1, stat = "identity", color="white") + facet_grid(~orig.ident)+ coord_polar("y", start=0) + scale_fill_manual(values = celltype.colv)

ggsave("pie_overlapped_chart.png", width = 6, height = 3, dpi = 300)
test <- subset(tnk, celltype == "CD8 Tem")
test$celltype %>% table
test <- subset(test, overlapped == "N")
test <- subset(test, orig.ident == "T92"| orig.ident == "T110")
test$celltype %>% table

Idents(test) <- "orig.ident"

test_marker <- FindMarkers(test, ident.1 = "T92", ident.2 = "T110", test.use = 'MAST')
test_marker <- test_marker %>% filter(p_val_adj < 0.05)

test <- ScaleData(test, features = rownames(test_marker))


test <- test@assays$RNA@scale.data %>% as.matrix()
library(pheatmap)

pheatmap(test,  cluster_cols = F, scale = "column",color=colorRampPalette(c("navy", "white", "red"))(50))

library(stringr)

colnames(test)<- str_split_fixed(test %>% colnames, "_", n = 2)[,1]

tnk_pair
tnk92 <- tnk_pair %>% subset(orig.ident == "T92")
tnk110 <- tnk_pair %>% subset(orig.ident == "T110")
tnk92$overlapped <- ifelse(tnk92$CTaa %in% tnk110$CTaa, "Overlapped", "Unique")
tnk110$overlapped <- ifelse(tnk110$CTaa %in% tnk92$CTaa, "Overlapped", "Unique")

tnk_pair <- merge(tnk92, tnk110)
tnk_pair <- subset(tnk_pair, celltype != "NK cells")

tnk_pair$overlapped %>% table
table(tnk_pair$orig.ident, tnk_pair$overlapped)

VlnPlot(tnk_pair, features = c("ZNF683", "ITGAE"), group.by = "overlapped")

tnk_overlapped <- subset(tnk_pair, overlapped == "Overlapped")
tnk_overlapped <- subset(tnk_overlapped, celltype != "NK cells")

VlnPlot(tnk_overlapped, features = c("ZNF683", "ITGAE"), group.by = "celltype" )

table(tnk_overlapped$celltype, tnk_overlapped$orig.ident)
tnk_overlapped@meta.data %>% View

meta <- meta[!is.na(meta$cloneType),]

