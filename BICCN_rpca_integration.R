rm(list = ls())
setwd("~/Desktop/Xinyu_Patchseq/")

library("Seurat")
library(readxl)
metadata <- as.data.frame(read_excel("reference_data/metadata.xlsx", 
                       sheet = "STable 1 - Whole Brain Metadata", 
                       col_types = c("text", "skip", "text",
                                     "skip", "text", "skip", "skip", "skip", 
                                     "skip", "text", "skip", "skip")))
metadata$cell.name = gsub("gw","GW",metadata$cell.name)
rownames(metadata)=metadata$cell.name
metadata = metadata[metadata$area=="PFC",-1]

library(stringr)
#metadata$cell = str_sub(metadata$cell.name,-16,-1)
#sum(metadata$cell%in%str_sub(colnames(reference),1,16))

data_dir = list.dirs(path = './reference_data', full.names = TRUE, recursive = TRUE)[-1]
seuobj.list = list()
for(i in 1:length(data_dir)){
  expression_matrix <- Read10X(data.dir = data_dir[i])
  cellname = intersect(
    paste(gsub("./reference_data/","",data_dir[i]),
          "_",str_sub(colnames(expression_matrix),1,16),sep=""),
    rownames(metadata))
  colnames(expression_matrix) = paste(gsub("./reference_data/","",data_dir[i]),
                                      "_",str_sub(colnames(expression_matrix),1,16),sep="")
  seuobj = CreateSeuratObject(counts = expression_matrix[,cellname])
  seuobj = AddMetaData(object = seuobj, metadata = metadata[cellname,3], col.name = 'celltype')
  seuobj.list[i] = seuobj
}
names(seuobj.list) = gsub("./reference_data/","",data_dir)
seuobj.list$GW22T_PFC = NULL
seuobj.list$mouse = NULL

seuobj.list <- lapply(X = seuobj.list, function(x) {
  x <- NormalizeData(x, verbose = TRUE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = seuobj.list)
seuobj.list <- lapply(X = seuobj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = TRUE)
})

anchors <- FindIntegrationAnchors(object.list = seuobj.list, anchor.features = features, verbose = TRUE, reduction = "rpca")
reference <- IntegrateData(anchorset = anchors, verbose = TRUE)
DefaultAssay(reference) <- "integrated"
reference <- ScaleData(reference, verbose = TRUE)

############################################

reference <- RunPCA(reference, verbose = TRUE)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:30, verbose = TRUE)
p1 <- DimPlot(reference, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(reference, reduction = "umap", group.by = "celltype")
p1 + p2

save(reference, file="BICCN_rpca_v2.rda")

FeaturePlot(reference,features = c("SOX2","EOMES","MKI67","HOPX","PDGFRA","AQP4","BCL11B","SATB2","NEUROD6","DLX6-AS1"))
##########
load("./BICCN_rpca_v2.rda")
human.query <- human_seuobj
human.anchors <- FindTransferAnchors(reference = reference, query = human.query,
                                        dims = 1:50, reference.reduction = "pca")
human.predictions <- TransferData(anchorset = human.anchors, refdata = reference$celltype,
                            dims = 1:50)
human.query <- AddMetaData(human.query, metadata = human.predictions)

##########
monkey.query <- monkey_seuobj
monkey.anchors <- FindTransferAnchors(reference = reference, query = monkey.query,
                                     dims = 1:30, reference.reduction = "pca")
monkey.predictions <- TransferData(anchorset = monkey.anchors, refdata = reference$celltype,
                                  dims = 1:30)
monkey.query <- AddMetaData(monkey.query, metadata = monkey.predictions)

##########

reference <- RunUMAP(reference, dims = 1:30, reduction = "pca", return.model = TRUE)
human.query <- MapQuery(anchorset = human.anchors, reference = reference, query = human.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

monkey.query <- MapQuery(anchorset = monkey.anchors, reference = reference, query = monkey.query,
                        refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")


p1 <- DimPlot(reference, reduction = "umap", group.by = "celltype") + ggtitle("Reference annotations (BICCN)")
p2 <- DimPlot(human.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transferred labels (Human)")
p3 <- DimPlot(monkey.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transferred labels (Monkey)")
p1+p2+p3



p4 <- DimPlot(human.query, reduction = "ref.umap", group.by = "pred2.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transferred labels (Human)")
p5 <- DimPlot(monkey.query, reduction = "ref.umap", group.by = "pred2.celltype", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transferred labels (Monkey)")
p1+p4+p5

FeaturePlot(human.query, reduction ="ref.umap",features = c("SOX2","EOMES","MKI67","HOPX","PDGFRA","AQP4","BCL11B","SATB2","NEUROD6","DLX6-AS1"))

save(human.query, file="human_query2.rda")
save(monkey.query, file="monkey_query2.rda")
