rm(list=ls())
setwd("~/Desktop/Xinyu_Patchseq/")
load("./shiny/data/BICCN_rpca.rda")
load("./shiny/data/human_query.rda")
load("./shiny/data/monkey_query.rda")

library(ggplot2)
library(gridExtra)

BICCN_umap = as.data.frame(reference@reductions[["umap"]]@cell.embeddings)
human_umap = as.data.frame(human.query@reductions[["ref.umap"]]@cell.embeddings)
monkey_umap = as.data.frame(monkey.query@reductions[["ref.umap"]]@cell.embeddings)

BICCN_umap$group = rep("BICCN",dim(BICCN_umap)[1])
BICCN_umap$celltype = factor(reference@meta.data[["celltype"]],
                             levels=c("Neuron","Endo","Astrocyte","Dividing","RG",
                                      "IPC","Interneuron","Vascular","Oligo",
                                      "Microglia","Outlier","NA"))

human_umap$group = rep("human",dim(human_umap)[1])
human_umap$celltype = factor(human.query@meta.data[["predicted.celltype"]],
                             levels=c("Neuron","Endo","Astrocyte","Dividing","RG",
                                      "IPC","Interneuron","Vascular","Oligo",
                                      "Microglia","Outlier","NA"))

monkey_umap$group = rep("monkey",dim(monkey_umap)[1])
monkey_umap$celltype = factor(monkey.query@meta.data[["predicted.celltype"]],
                              levels=c("Neuron","Endo","Astrocyte","Dividing","RG",
                                       "IPC","Interneuron","Vascular","Oligo",
                                       "Microglia","Outlier","NA"))

colnames(human_umap) = colnames(BICCN_umap)
colnames(monkey_umap) = colnames(BICCN_umap)
umaps = rbind(BICCN_umap,human_umap,monkey_umap)

umaps$human_celltype = factor(c(rep("NA",81337),as.character(human_umap$celltype),rep("NA",270)),
                              levels=c("Neuron","Endo","Astrocyte","Dividing","RG",
                                       "IPC","Interneuron","Vascular","Oligo",
                                       "Microglia","Outlier","NA"))
umaps$monkey_celltype = factor(c(rep("NA",81337),rep("NA",145),as.character(monkey_umap$celltype)),
                               levels=c("Neuron","Endo","Astrocyte","Dividing","RG",
                                        "IPC","Interneuron","Vascular","Oligo",
                                        "Microglia","Outlier","NA"))

col = c("#A83337","#FF924C","#FFCA3A","#C5CA30","#8AC926","#52A675","#5DAAD9","#4267AC","#6A4C93","#7D7F7F","#E5E8E8")

col2 = c("#A83337","#FF924C","#FFCA3A","#C5CA30","#8AC926","#52A675","#5DAAD9","#4267AC","#6A4C93","#43256B","#E5E8E8","#E5E8E8")

################################# fig version1 #################################
p0 = ggplot(umaps, aes(x=UMAP_1, y=UMAP_2,color=celltype)) +
  geom_point(size = 2) + 
  scale_color_manual(values = col)+
  #geom_point(aes(colour = `Days in vitro (DIV)`))+
  #scale_colour_gradient(low = "grey", high = "blue") + 
  ggtitle("BICCN (Reference)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


p1 = ggplot(umaps, aes(x=UMAP_1, y=UMAP_2,color=celltype)) +
  geom_point(size = 0.15, alpha = 0.4) + 
  scale_color_manual(values = col)+
  #geom_point(aes(colour = `Days in vitro (DIV)`))+
  #scale_colour_gradient(low = "grey", high = "blue") + 
  ggtitle("BICCN (Reference)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),legend.position = "none")

p2 = ggplot(umaps[umaps$group!="monkey",], aes(x=UMAP_1, y=UMAP_2,
                                         color = human_celltype, 
                                         shape = group, 
                                         size = group)) +
  geom_point(alpha = 0.85) + 
  scale_color_manual(values = c(col2[c(1,3,4,5,7,9,10)],"#E5E8E8"))+
  scale_shape_manual(values = c(20, 15)) +
  scale_size_manual(values = c(0.15, 1.2)) +
  #geom_point(aes(colour = `Days in vitro (DIV)`))+
  #scale_colour_gradient(low = "grey", high = "blue") + 
  ggtitle("Human") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+ guides(shape = "none",color = "none", size = "none")

p3=ggplot(umaps[umaps$group!="human",], aes(x=UMAP_1, y=UMAP_2,
                                         color = monkey_celltype, 
                                         shape = group, 
                                         size = group)) +
  geom_point(alpha = 0.85) + 
  scale_color_manual(values = col2[c(1,3,4,5,9,10,12)])+
  scale_shape_manual(values = c(20, 15)) +
  scale_size_manual(values = c(0.15, 1.2)) +
  #geom_point(aes(colour = `Days in vitro (DIV)`))+
  #scale_colour_gradient(low = "grey", high = "blue") + 
  ggtitle("Monkey") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+ guides(shape = "none",color = "none", size = "none")
grid.arrange(p1, p2, p3, nrow = 1)

## save as .tiff
tiff("annotation1.tiff", units="in", width=15, height=5, res=600)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

## save as .pdf
pdf("annotation.pdf", width=15, height=5)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

################################# fig version2 #################################
cells

umaps2=umaps[umaps$group=="BICCN"|rownames(umaps)%in%cells,]

p4 = ggplot(umaps2[umaps2$group!="monkey",], aes(x=UMAP_1, y=UMAP_2,
                                               color = human_celltype, 
                                               shape = group, 
                                               size = group)) +
  geom_point(alpha = 0.85) + 
  scale_color_manual(values = c(col2[c(1,3,5,9,10)],"#E5E8E8"))+
  scale_shape_manual(values = c(20, 15)) +
  scale_size_manual(values = c(0.15, 1.2)) +
  #geom_point(aes(colour = `Days in vitro (DIV)`))+
  #scale_colour_gradient(low = "grey", high = "blue") + 
  ggtitle("Human") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+ guides(shape = "none",color = "none", size = "none")

p5 = ggplot(umaps2[umaps2$group!="human",], aes(x=UMAP_1, y=UMAP_2,
                                            color = monkey_celltype, 
                                            shape = group, 
                                            size = group)) +
  geom_point(alpha = 0.85) + 
  scale_color_manual(values = col2[c(1,3,4,5,9,10,12)])+
  scale_shape_manual(values = c(20, 17)) +
  scale_size_manual(values = c(0.15, 1.2)) +
  #geom_point(aes(colour = `Days in vitro (DIV)`))+
  #scale_colour_gradient(low = "grey", high = "blue") + 
  ggtitle("Monkey") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+ guides(shape = "none",color = "none", size = "none")
grid.arrange(p1, p4, p5, nrow = 1)


tiff("annotation2.tiff", units="in", width=15, height=5, res=600)
grid.arrange(p1, p4, p5, nrow = 1)
dev.off()

tiff("annotation0.tiff", units="in", width=7, height=5, res=600)
grid.arrange(p0, nrow = 1)
dev.off()

################################# fig version3 #################################

celltype_all = umaps2$celltype
celltype_all[umaps2$group=="BICCN"]="NA"
umaps2$celltype_all = celltype_all

p6 = ggplot(umaps2, aes(x=UMAP_1, y=UMAP_2,
                        color = celltype_all, 
                        shape = group, 
                        size = group)) +
  geom_point(alpha = 0.85) + 
  scale_color_manual(values = c(col2[c(1,3,4,5,9,10)],"#E5E8E8"))+
  scale_shape_manual(values = c(20, 15, 17)) +
  scale_size_manual(values = c(0.15, 1.2, 1.2)) +
  #geom_point(aes(colour = `Days in vitro (DIV)`))+
  #scale_colour_gradient(low = "grey", high = "blue") + 
  ggtitle("") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))#+ guides(shape = "none",color = "none", size = "none")

grid.arrange(p1, p6, nrow = 1)


tiff("annotation3.tiff", units="in", width=10, height=5, res=600)
grid.arrange(p1, p6, nrow = 1)
dev.off()

## eFeature plot

######################################################################################################
e.human$PCD = paste("PCD ",e.human$`Gestation Day (GD)`,".",e.human$`Days in vitro (DIV)`,sep="")
unique(e.human$PCD)
#e.human$PCD = factor(e.human$PCD,levels = c("PCD 132.0","PCD 132.7","PCD 132.14","PCD 132.18","PCD 132.19","PCD 132.20","PCD 132.21"))
e.human = e.human[colnames(human.query),]
human_umap$PCD = e.human$PCD




######################################################################################################
e.monkey$PCD = paste("PCD ",e.monkey$`Gestation Day (GD)`,".",e.monkey$`Days in vitro (DIV)`,sep="")
#unique(e.monkey$PCD)
#e.monkey$PCD = factor(e.monkey$PCD,levels = c("PCD 90.0","PCD 90.6","PCD 90.10","PCD 90.16","PCD 90.22", 
#                                              "PCD 96.0","PCD 96.16",
#                                              "PCD 100.0","PCD 100.12","PCD 100.19",
#                                              "PCD 106.0",
#                                              "PCD 112.0",
#                                              "PCD 155.0"))
e.monkey = e.monkey[colnames(monkey.query),]
monkey_umap$PCD = e.monkey$PCD

######
umaps$PCD = c(rep("NA",81337),human_umap$PCD,monkey_umap$PCD)
umaps$human_PCD = c(rep("NA",81337),human_umap$PCD,rep("NA",270))
umaps$human_PCD = factor(umaps$human_PCD, levels = c("NA","PCD 132.0","PCD 132.7","PCD 132.14","PCD 132.18","PCD 132.19","PCD 132.20","PCD 132.21"))

umaps$monkey_PCD = c(rep("NA",81337),rep("NA",145),monkey_umap$PCD)
umaps$monkey_PCD = factor(umaps$monkey_PCD, levels = c("NA","PCD 90.0","PCD 90.6","PCD 90.10","PCD 90.16","PCD 90.22", 
                                                     "PCD 96.0","PCD 96.16",
                                                     "PCD 100.0","PCD 100.12","PCD 100.19",
                                                     "PCD 106.0",
                                                     "PCD 112.0",
                                                     "PCD 155.0"))


p4 = ggplot(umaps[umaps$group!="monkey",], aes(x=UMAP_1, y=UMAP_2,
                       color=human_PCD,size = group)) +
  geom_point(alpha = 0.85) + 
  scale_color_manual(values = c('#E5E8E8',rev(brewer.pal(7, "Greens"))))+
  scale_size_manual(values = c(0.15, 1.2))+
  labs(title="PCD_Human") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(),plot.title = element_text(hjust = 0.5))+ guides( size = "none")

p5 = ggplot(umaps[umaps$group!="human",], aes(x=UMAP_1, y=UMAP_2,
                                          color=monkey_PCD,size = group)) +
  geom_point(alpha = 0.85) + 
  scale_color_manual(values = c('#E5E8E8',
                                rev(brewer.pal(5, "Reds")),
                                "#E7E40D","#EFED91",
                                rev(brewer.pal(3, "Greens")),
                                "#3182BD",
                                "#756BB1",
                                "#523650"))+
  scale_size_manual(values = c(0.15, 1.2))+
  labs(title="PCD_Monkey") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(),plot.title = element_text(hjust = 0.5))+ guides( size = "none")

tiff("pcd.tiff", units="in", width=12, height=5, res=300)
grid.arrange(p4, p5, nrow = 1)
dev.off()

tiff("pcd1.tiff", units="in", width=6, height=5, res=300)
p4
dev.off()

tiff("pcd2.tiff", units="in", width=6, height=5, res=300)
p5
dev.off()

(table(human_umap$celltype,human_umap$PCD)["Neuron",]/table(human_umap$PCD))*100

(table(monkey_umap$celltype,monkey_umap$PCD)["Neuron",]/table(monkey_umap$PCD))*100

genelist

tiff("feature_human.tiff", units="in", width=15, height=43.5, res=300)
FeaturePlot(human.query,features=genelist)
dev.off()

tiff("feature_monkey.tiff", units="in", width=15, height=43.5, res=300)
FeaturePlot(monkey.query,features=genelist)
dev.off()
