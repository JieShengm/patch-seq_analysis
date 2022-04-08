## Trajectory
library(monocle)
library(dplyr)
library(Seurat)
library(RColorBrewer)

rm(list=ls())
setwd("~/Desktop/Xinyu_Patchseq/")
load("./neuron_monkey.rda")

### 1 Monkey neurons
### 2 Monkey neurons (Nutrite == "Yes")
#### 2.2 For non-cultured cells 
# 1 load filtered dataset from seurat object
seu.0=seu[,rownames(e.monkey)[e.monkey$`Neurite?`=="Yes" &e.monkey$`Days in vitro (DIV)`!= 0]]
seu.0=NormalizeData(seu.0)
seu.0=FindVariableFeatures(seu.0,nfeatures = 2000)
seu.0=ScaleData(seu.0)
seu.0 <- RunPCA(seu.0, npcs = 30, verbose = FALSE)
seu.0 <- RunUMAP(seu.0, reduction = "pca", dims = 1:30)

mat.0 = as.matrix(seu.0@assays$RNA@counts[VariableFeatures(seu.0),])
barcodes = colnames(seu.0)
features = data.frame(row.names=VariableFeatures(seu.0))
group = e.monkey[barcodes,]$PCD
gestation = factor(e.monkey[barcodes,]$`Gestation Day (GD)`)
meta.0 = data.frame(barcodes, group, gestation)
pd <- new("AnnotatedDataFrame", data = meta.0)
fd <- new("AnnotatedDataFrame", data = features)
fd$gene_short_name = VariableFeatures(seu.0)
colnames(mat.0) = rownames(barcodes)
rownames(mat.0) = rownames(features)
## create monocle object
mono_obj.0 <- newCellDataSet(as(mat.0,  "sparseMatrix"),
                             phenoData = pd, 
                             featureData = fd,
                             expressionFamily=negbinomial())

# calculate basic summary statistics of the expressino data
mono_obj.0 <- estimateSizeFactors(mono_obj.0)
mono_obj.0 <- estimateDispersions(mono_obj.0)

mono_obj.0 <- detectGenes(mono_obj.0, min_expr = 1)

diff_test_res.0 <- differentialGeneTest(mono_obj.0, fullModelFormulaStr = "~group")
ordering_genes <- row.names(subset(diff_test_res.0, qval < 0.05))
mono_obj.0 <- setOrderingFilter(mono_obj.0, ordering_genes)
mono_obj.0 <- reduceDimension(mono_obj.0, max_components = 3, method = 'DDRTree')
mono_obj.0 <- orderCells(mono_obj.0)
#gestation_col = c("#CB181D","#E7E40D","#31A354","#3182BD","#756BB1","#523650")
plot_cell_trajectory(mono_obj.0, color_by = "group",cell_size = 3,show_tree = TRUE,show_branch_points = FALSE)# + scale_color_manual(values=gestation_col)

######################

# 1 load filtered dataset from seurat object
seu.1=seu[,rownames(e.monkey)[e.monkey$`Neurite?`=="Yes" &e.monkey$`Gestation Day (GD)`==90]]
seu.1=NormalizeData(seu.1)
seu.1=FindVariableFeatures(seu.1)
seu.1=ScaleData(seu.1)
seu.1 <- RunPCA(seu.1, npcs = 30, verbose = FALSE)
seu.1 <- RunUMAP(seu.1, reduction = "pca", dims = 1:30)

mat.1 = as.matrix(seu.1@assays$RNA@counts)
barcodes = colnames(seu.1)
features = data.frame(row.names=rownames(seu.1))
group = e.monkey[barcodes,]$PCD
day = factor(e.monkey[barcodes,]$`Days in vitro (DIV)`)
meta.1 = data.frame(barcodes, group, day)
pd <- new("AnnotatedDataFrame", data = meta.1)
fd <- new("AnnotatedDataFrame", data = features)
fd$gene_short_name = rownames(seu.1)
colnames(mat.1) = rownames(barcodes)
rownames(mat.1) = rownames(features)
## create monocle object
mono_obj.1 <- newCellDataSet(as(mat.1,  "sparseMatrix"),
                             phenoData = pd, 
                             featureData = fd,
                             expressionFamily=negbinomial())

# calculate basic summary statistics of the expressino data
mono_obj.1 <- estimateSizeFactors(mono_obj.1)
mono_obj.1 <- estimateDispersions(mono_obj.1)

mono_obj.1 <- detectGenes(mono_obj.1, min_expr = 1)

diff_test_res <- differentialGeneTest(mono_obj.1, fullModelFormulaStr = "~day")
ordering_genes1 <- row.names(subset(diff_test_res, qval < 0.2))
mono_obj.1 <- setOrderingFilter(mono_obj.1, ordering_genes1)

mono_obj.1 <- reduceDimension(mono_obj.1, max_components = 2, method = 'DDRTree')
mono_obj.1 <- orderCells(mono_obj.1)

plot_cell_trajectory(mono_obj.1, color_by = "day",show_tree = FALSE,show_branch_points = FALSE) + scale_color_manual(values=rev(brewer.pal(4, "Reds")))

