library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
olbre.path <- paste0(wd.path, "olbrecht/")
dir.create(olbre.path)
setwd(olbre.path)

panc <- readRDS("d:/project_directory_name/ovary_project/Seurat/final/olbrecht/panc_tnk.rds")

panc$orig.ident %>% table
Idents(panc) <- "orig.ident"
#"BT1304", "BT1305",
olbrecht <- subset(panc, idents = c("BT1303", "BT1306", "BT1307"))
olbrecht$tumor %>% table

olbemarkers <- list() 
olbemarkers$activated <- c("CXCL13", "CD200", "ICOS")
olbemarkers$resting <- c("CCR7", "IL7R", "GPR183", "LMNA")
olbemarkers$regulatory <- c("IL2RA","TNFRSF9", "FOXP3", "CTLA4")


olbrecht<- AddModuleScore(olbrecht, features = olbemarkers, name = c("activated", "resting", "regulatory"))
Idents(olbrecht) <- "tumor"

extuolb <- subset(olbrecht, idents = "1")

olbrecht@meta.data %>% colnames

olmeta<- extuolb@meta.data
olmeta$TumorSite %>% table
olmeta$TumorSite <- factor(olmeta$TumorSite, levels = c("Ovarium", "Omentum", "Peritoneum"))
TumorSite.color <- c(brewer.pal(3, "Set2"))
olmeta <- olmeta %>% arrange(desc(TumorSite))

p1 <- ggplot(olmeta, aes(x = TumorSite, y = resting2, fill = TumorSite)) + 
  geom_boxplot(outlier.color = "black", outlier.size = 0.2, size= 0.2) + theme_classic() + theme(axis.text.x = element_blank()) + labs(x = " ", y = "Resting_score") + NoLegend() + scale_fill_manual(values = TumorSite.color)

p2 <- ggplot(olmeta, aes(x = TumorSite, y = activated1, fill = TumorSite)) + geom_boxplot(outlier.color = "black", outlier.size = 0.2, size= 0.2 ) + theme_classic() +
  theme(axis.text.x = element_blank()) + labs(x = " ",y = "Activated_score") + NoLegend() + scale_fill_manual(values = TumorSite.color)

p3 <- ggplot(olmeta, aes(x = TumorSite, y = regulatory3, fill = TumorSite)) + geom_boxplot(outlier.color = "black", outlier.size = 0.2, size= 0.2) + theme_classic() + 
  theme(axis.text.x = element_blank()) + labs(x = " ", y ="Regulatory_score") + scale_fill_manual(values = TumorSite.color)

p1 + p2 + p3 + patchwork::plot_layout(ncol =  3)

ggsave("olbrecht_et_al.png", width = 6, height = 2.5, dpi = 300)
ggsave("olbrecht_et_al.pdf", width = 6, height = 2.5, dpi = 600)
