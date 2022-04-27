rm(list=ls())
load("./neuron_monkey.rda")

cells=rownames(e.monkey)[e.monkey$`Days in vitro (DIV)`==0&e.monkey$`Neurite?`=="Yes"]

seu = seu[,cells]
seu = NormalizeData(seu)
seu = ScaleData(seu,features = rownames(seu))

mat = as.matrix(seu[,cells]@assays$RNA@data)
colnames(e.monkey[cells,c(8:10,13)])

gene.ind_spearman = apply((cor(t(mat),e.monkey[cells,c(8:10,13)],method="spearman")),2,function(x){order(abs(x),decreasing=TRUE)[1:30]})
corr_spearman=apply((cor(t(mat),e.monkey[cells,c(8:10,13)],method="spearman")),2,function(x){sort(abs(x),decreasing=TRUE)[1:30]})
corr_gene_spearman = apply(gene.ind_spearman,2,function(x){rownames(seu[,cells])[x]})

gene.ind_pearson = apply((cor(t(mat),e.monkey[cells,c(8:10,13)])),2,function(x){order(abs(x),decreasing=TRUE)[1:30]})
corr_pearson=apply((cor(t(mat),e.monkey[cells,c(8:10,13)])),2,function(x){sort(abs(x),decreasing=TRUE)[1:30]})
corr_gene_pearson = apply(gene.ind_pearson,2,function(x){rownames(seu[,cells])[x]})


cor(t(mat[corr_gene_spearman[,1],]),e.monkey[cells,8],method="spearman")
cor(t(mat[corr_gene_spearman[,2],]),e.monkey[cells,9],method="spearman")
cor(t(mat[corr_gene_spearman[,3],]),e.monkey[cells,10],method="spearman")
cor(t(mat[corr_gene_spearman[,4],]),e.monkey[cells,13],method="spearman")
#####################
write.csv(data.frame(corr_gene_spearman[,1],correlation_1 = cor(t(mat[corr_gene_spearman[,1],]),e.monkey[cells,8],method="spearman"),
                     corr_gene_spearman[,2],correlation_1 = cor(t(mat[corr_gene_spearman[,2],]),e.monkey[cells,9],method="spearman"),
                     corr_gene_spearman[,3],correlation_1 = cor(t(mat[corr_gene_spearman[,3],]),e.monkey[cells,10],method="spearman"),
                     corr_gene_spearman[,4],correlation_1 = cor(t(mat[corr_gene_spearman[,4],]),e.monkey[cells,13],method="spearman"),row.names = 1:30),
          file="correlated_genes(spearman_correlation).csv")

library(Seurat)
library(patchwork)
library(ggplot2)

# "Membrane Capacitance (pF)"       
# "Input Resistance (M\u03a9)"     
# "Resting Membrane Potential (mV)" 
# "INa(pA/pF)"           

VlnPlot(seu[,cells], corr_gene_spearman[,1], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + 
  ggtitle("Top 30 genes correlated with Membrane Capacitance")

VlnPlot(seu[,cells], corr_gene_spearman[,2], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + ggtitle("Top 30 genes correlated with Input Resistance")

VlnPlot(seu[,cells], corr_gene_spearman[,3], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + ggtitle("Top 30 genes correlated with Resting Membrane Potential")

VlnPlot(seu[,cells], corr_gene_spearman[,4], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + ggtitle("Top 30 genes correlated with INa")


VlnPlot(seu[,cells], corr_gene_spearman[,3], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + ggtitle("Top 30 genes correlated with Resting Membrane Potential")

DotPlot(seu[,cells], features = corr_gene_spearman[,1]) + RotatedAxis()+coord_flip()
DoHeatmap(seu[,cells], features = corr_gene_spearman[,4])
markers = FindAllMarkers(seu)
############################
write.csv(data.frame(corr_gene_pearson[,1],correlation_1 = cor(t(mat[corr_gene_pearson[,1],]),e.monkey[cells,8]),
                     corr_gene_pearson[,2],correlation_1 = cor(t(mat[corr_gene_pearson[,2],]),e.monkey[cells,9]),
                     corr_gene_pearson[,3],correlation_1 = cor(t(mat[corr_gene_pearson[,3],]),e.monkey[cells,10]),
                     corr_gene_pearson[,4],correlation_1 = cor(t(mat[corr_gene_pearson[,4],]),e.monkey[cells,13]),row.names = 1:30),
          file="correlated_genes(pearson_correlation).csv")

VlnPlot(seu[,cells], corr_gene_pearson[,1], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + 
  ggtitle("Top 30 genes correlated with Membrane Capacitance")

VlnPlot(seu[,cells], corr_gene_pearson[,2], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + ggtitle("Top 30 genes correlated with Input Resistance")

VlnPlot(seu[,cells], corr_gene_pearson[,3], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + ggtitle("Top 30 genes correlated with Resting Membrane Potential")

VlnPlot(seu[,cells], corr_gene_pearson[,4], stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = rel(1),face = "bold", angle = 0)) + ggtitle("Top 30 genes correlated with INa")
