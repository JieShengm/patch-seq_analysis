rm(list = ls())
setwd("~/Desktop/Xinyu_Patchseq/")

library("Seurat")
library(RColorBrewer)
metadata <- as.data.frame(read.csv("reference_data/mouse/metaData_scDevSC.csv",sep='\t',row.names = 1))
cells = rownames(metadata)[metadata$biosample_id=="P4"]

data_dir <- 'reference_data/mouse/'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
expression_matrix=expression_matrix[,colnames(expression_matrix)%in%cells]
reference = CreateSeuratObject(counts = expression_matrix)
reference = AddMetaData(object = reference, metadata = metadata[cells,], col.name = colnames(metadata[cells,]))

expression_matrix = read.table("./reference_data/mouse/wt36ko36.txt")
query = CreateSeuratObject(counts = expression_matrix)
  
reference <- NormalizeData(reference, verbose = TRUE)
reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
reference = ScaleData(reference)

query <- NormalizeData(query , verbose = TRUE)
query  <- FindVariableFeatures(query , selection.method = "vst", nfeatures = 2000, verbose = TRUE)
query  = ScaleData(query )

reference <- RunPCA(reference, verbose = TRUE)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:50, verbose = TRUE)
p1 = DimPlot(reference, reduction = "umap", group.by = "New_cellType") 
p1

anchors0 <- FindTransferAnchors(reference = reference, query = query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors0, refdata = reference$New_cellType,
                            dims = 1:30)
query <- AddMetaData(query, metadata = predictions)

reference <- RunUMAP(reference, dims = 1:50, reduction = "pca", return.model = TRUE)
query <- MapQuery(anchorset = anchors0, reference = reference, query = query,
                           refdata = list(celltype = "New_cellType"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(reference, reduction = "umap", group.by = "New_cellType")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) 
p1+p2

table(query$predicted.celltype)



################################################
## 249 low quality cells
reference2 = reference[,reference$New_cellType!="Low quality cells"]
reference2 <- NormalizeData(reference2, verbose = TRUE)
reference2 <- FindVariableFeatures(reference2, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
reference2 = ScaleData(reference2)


reference2 <- RunPCA(reference2, verbose = TRUE)
reference2 <- RunUMAP(reference2, reduction = "pca", dims = 1:50, verbose = TRUE)
p1 = DimPlot(reference2, reduction = "umap", group.by = "New_cellType") 
p1

anchors2 <- FindTransferAnchors(reference = reference2, query = query2,
                                dims = 1:30, reference.reduction = "pca")
predictions2 <- TransferData(anchorset = anchors2, refdata = reference2$New_cellType,
                            dims = 1:30,k.weight = 30)
query2 <- AddMetaData(query2, metadata = predictions)

reference2 <- RunUMAP(reference2, dims = 1:30, reduction = "pca", return.model = TRUE)
query2 <- MapQuery(anchorset = anchors2, reference = reference2, query = query2,
                  refdata = list(celltype = "New_cellType"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(reference2, reduction = "umap", group.by = "New_cellType")
p2 <- DimPlot(query2, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 3, repel = TRUE) 
p1+p2

table(query$predicted.celltype)
