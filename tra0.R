rm(list=ls())

setwd("~/Desktop/Xinyu_Patchseq/")
load("./shiny/data/monkey_query.rda")

library(slingshot)
library(Seurat)
################################################################################
seu = (monkey.query[,monkey.query@meta.data[["predicted.id"]]=="Neuron"])
seu = AddMetaData(seu, metadata = e.monkey$PCD, col.name = "PCD")
sce = as.SingleCellExperiment(seu)
table(seu@meta.data$PCD)

sce <- runPCA(sce, ncomponents = 50)
sce$PC1 = reducedDim(sce)[,1]
sce$PC2 = reducedDim(sce)[,2]
  
ggplot(as.data.frame(colData(sce)), aes(x = PC1, y = PC2, color = PCD)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = pcd_col) + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

sce$pseudotime_PC1 <- rank(sce$PC1)  # rank cells by their PC1 score
ggplot(as.data.frame(colData(sce)), aes(x = pseudotime_PC1, y = PCD, 
                                             colour = PCD)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = pcd_col)+
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

################################################################################
library("destiny")
des_sce <- logcounts(sce)  # access log-transformed counts matrix
cellLabels <- sce$PCD
colnames(des_sce) <- cellLabels

# Make a diffusion map.
# Optional: Try different sigma values when making diffusion map.
dm <- DiffusionMap(as.matrix(t(des_sce)), sigma = "local")  # use local option to set sigma
sigmas <- find_sigmas(as.matrix(t(des_sce)), verbose = TRUE)  # find optimal sigma
dm <- DiffusionMap(as.matrix(t(des_sce)), sigma = optimal_sigma(sigmas))  

# Plot diffusion component 1 vs diffusion component 2 (DC1 vs DC2). 
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = sce$PCD,
                  gestation = gestation,
                  culture = culture)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_manual(values = pcd_col) + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_manual(values = c(rev(brewer.pal(5, "Reds")),
                                               "#E7E40D","#EFED91",
                                               rev(brewer.pal(3, "Greens")),
                                               "#3182BD",
                                               "#756BB1",
                                               "#523650")) + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

# Try plotting higher diffusion components against one another.

# Next, let us use the first diffusion component (DC1) as a measure of pseudotime.
# How does the separation by cell stage look?
sce$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])    # rank cells by their dpt
ggplot(as.data.frame(colData(sce)), 
       aes(x = pseudotime_diffusionmap, 
           y = PCD, colour = PCD)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = pcd_col) + theme_classic() +
  xlab("Diffusion component 1 (DC1)") + ylab("Timepoint") +
  ggtitle("Cells ordered by DC1")

sce$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])    # rank cells by their dpt
ggplot(as.data.frame(colData(sce)), 
       aes(x = pseudotime_diffusionmap, 
           y = gestation, colour = gestation)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = pcd_col) + theme_classic() +
  xlab("Diffusion component 1 (DC1)") + ylab("Timepoint") +
  ggtitle("Cells ordered by DC1")

sce$pseudotime_diffusionmap2 <- rank(eigenvectors(dm)[,2])    # rank cells by their dpt
ggplot(as.data.frame(colData(sce)), 
       aes(x = pseudotime_diffusionmap2, 
           y = PCD, colour = PCD)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = pcd_col) + theme_classic() +
  xlab("Diffusion component 1 (DC2)") + ylab("Timepoint") +
  ggtitle("Cells ordered by DC2")

plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')

#####
pca = reducedDim(sce,"PCA")
# What happens if you run the diffusion map on the PCs? Why would one do this?
rownames(pca) <- cellLabels
dm <- DiffusionMap(pca)

# Diffusion pseudotime calculation. 
# Set index or tip of pseudotime calculation to be a zygotic cell (cell 268). 
dpt <- DPT(dm, tips = 1)

# Plot DC1 vs DC2 and color the cells by their inferred diffusion pseudotime.
# We can accesss diffusion pseudotime via dpt$dpt.
df <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2], 
                 dptval = dpt$dpt, cell_type2 = sce$PCD)
p1 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = dptval))
p2 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = cell_type2))
p <- plot_grid(p1, p2)
p

# Plot diffusion pseudotime vs timepoint. 
# Which separates the data better, DC1 or diffusion pseudotime?
sce$pseudotime_dpt <- rank(dpt$dpt) 
ggplot(as.data.frame(colData(sce)), 
       aes(x = pseudotime_dpt, 
           y = PCD, colour = PCD)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = pcd_col) + theme_classic() +
  xlab("Diffusion map pseudotime (dpt)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")

################################################################################
library(slingshot)
library(Seurat)

# Read the Slingshot documentation (?slingshot) and then run Slingshot below. 
# Given your understanding of the algorithm and the documentation, what is one 
# major set of parameters we omitted here when running Slingshot?
sce <- slingshot(sce, reducedDim = 'PCA')  # no clusters

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

# Plot Slingshot pseudotime vs cell stage. 
ggplot(as.data.frame(colData(sce)), aes(x = sce$slingPseudotime_1, y = PCD, 
                                             colour = PCD)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_manual(values = pcd_col) + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# Cluster cells using the Seurat workflow below.
gcdata <- CreateSeuratObject(counts = counts(sce), min.cells = 0, min.genes = 0, project = "slingshot")
gcdata <- NormalizeData(object = gcdata)
gcdata <- FindVariableFeatures(object = gcdata, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.5)
gcdata <- ScaleData(object = gcdata, do.center = T, do.scale = F)
gcdata <- RunPCA(object = gcdata, pc.genes = VariableFeatures(gcdata), do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5)
gcdata <- FindClusters(object = gcdata, reduction = "pca", dims = 20, 
                       resolution = 0.6)

# Add clustering information from Seurat to the deng_SCE object
# Then run Slingshot using these cluster assignments.
#sce$slingPseudotime_seu <- NULL  # remove old slingshot pseudotime data
colData(sce)$Seurat_clusters <- as.character(gcdata@ident)  # go from factor to character
deng_SCE <- slingshot(deng_SCE, clusterLabels = 'Seurat_clusters', reducedDim = 'PCA')

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(deng_SCE)$PCA, col = colors[cut(deng_SCE$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(deng_SCE), lwd=2)

# Plot Slingshot pseudotime vs cell stage. 
ggplot(as.data.frame(colData(deng_SCE)), aes(x = slingPseudotime_1, y = cell_type2, 
                                             colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
ggsave(paste0(mydir, "/pseudotime_slingshot.png"))

################################################################################

library(monocle)
library(dplyr)

# 1 load filtered dataset from seurat object
cell_remained = colnames(seu)
## extract clusters detected by seurat
timepoint = as.data.frame(seu@meta.data[["PCD"]],row.names = cell_remained)
colnames(timepoint) = "timepoint"

# 2 create monocle object

seu.0=seu[,e.monkey[barcodes,]$`Days in vitro (DIV)`==0]
seu.0=NormalizeData(seu.0)
seu.0=FindVariableFeatures(seu.0)
seu.0=ScaleData(seu.0)


mat = as.matrix(seu.0@assays$RNA@counts[VariableFeatures(seu.0),])
barcodes = colnames(seu.0)
features = data.frame(row.names=rownames(seu.0))
group = e.monkey[barcodes,]$PCD
gestation = factor(e.monkey[barcodes,]$`Gestation Day (GD)`)

meta = data.frame(barcodes, group, gestation)
pd <- new("AnnotatedDataFrame", data = meta)
fd <- new("AnnotatedDataFrame", data = features)
fd$gene_short_name = rownames(seu)
colnames(mat) = rownames(barcodes)
rownames(mat) = rownames(features)
## create monocle object
mono_obj <- newCellDataSet(as(mat,  "sparseMatrix"),
                           phenoData = pd, 
                           featureData = fd,
                           expressionFamily=negbinomial())

# calculate basic summary statistics of the expressino data
mono_obj <- estimateSizeFactors(mono_obj)
mono_obj <- estimateDispersions(mono_obj)

mono_obj<- detectGenes(mono_obj, min_expr = 1)
expressed_genes <- row.names(subset(fData(mono_obj), num_cells_expressed >= 10))

############
diff_test_res <- differentialGeneTest(mono_obj[expressed_genes,], fullModelFormulaStr = "~timepoint")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
mono_obj <- setOrderingFilter(mono_obj, ordering_genes)

mono_obj <- reduceDimension(mono_obj, max_components = 2, method = 'DDRTree')
mono_obj <- orderCells(mono_obj)

plot_cell_trajectory(mono_obj, color_by = "timepoint")
plot_cell_trajectory(mono_obj, color_by = "timepoint") +
  facet_wrap(~State, nrow = 1)

plot_cell_trajectory(mono_obj, color_by = "gestation")
plot_cell_trajectory(mono_obj, color_by = "gestation") +
  facet_wrap(~State, nrow = 1)

plot_cell_trajectory(mono_obj, color_by = "culture")
plot_cell_trajectory(mono_obj, color_by = "culture") +
  facet_wrap(~State, nrow = 1)

plot_cell_trajectory(mono_obj, color_by = "Pseudotime")
plot_cell_trajectory(mono_obj, color_by = "Pseudotime") +
  facet_wrap(~State, nrow = 1)


group_vec <- levels(gestation)
group_cols <- brewer.pal(7, "Reds")
names(group_cols) <- group_vec
plot_complex_cell_trajectory(mono_obj[, ], color_by = 'gestation', show_branch_points = T, cell_size = 1, cell_link_size = 0.3) + facet_wrap(~culture, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = group_cols, name = "gestation")

group_vec <- levels(culture)
group_cols <- brewer.pal(7, "Spectral")
names(group_cols) <- group_vec
plot_complex_cell_trajectory(mono_obj[, ], color_by = 'culture', show_branch_points = T, cell_size = 1, cell_link_size = 0.3) + facet_wrap(~gestation, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = group_cols, name = "culture")

################################################################################
importCDS(data_to_be_imported)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")