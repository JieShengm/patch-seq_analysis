
load("./shiny/data/monkey_query.rda")
seu = (monkey.query[,monkey.query@meta.data[["predicted.id"]]=="Neuron"])

e.monkey = e.monkey[colnames(seu),]
seu = NormalizeData(seu)
seu = FindVariableFeatures(seu)
seu = ScaleData(seu)
seu = seu[VariableFeatures(seu),]
Idents(seu) = gestation
markers=FindMarkers(seu,ident.1 = "90",ident.2 = "155")
markers = rownames(markers)[markers$p_val_adj<0.1]


mat = as.matrix(seu[markers,]@assays$RNA@counts)
barcodes = colnames(seu)
#features = data.frame(row.names=rownames(seu))
features = data.frame(row.names = markers)
gestation = factor(e.monkey[barcodes,]$`Gestation Day (GD)`)
culture = factor(e.monkey[barcodes,]$`Days in vitro (DIV)`)

gestation_num = e.monkey[barcodes,]$`Gestation Day (GD)`

meta = data.frame(barcodes, gestation,culture)
pd <- new("AnnotatedDataFrame", data = meta)
fd <- new("AnnotatedDataFrame", data = features)
#fd$gene_short_name = rownames(seu)
fd$gene_short_name = markers
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
diff_test_res <- differentialGeneTest(mono_obj[expressed_genes,], fullModelFormulaStr = "~gestation")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
mono_obj <- setOrderingFilter(mono_obj, ordering_genes)

mono_obj <- reduceDimension(mono_obj, max_components = 2, method = 'DDRTree')
mono_obj <- orderCells(mono_obj)

g90=mono_obj[,gestation_num==90]

plot_cell_trajectory(mono_obj[,gestation_num==90], color_by = "gestation")+
  scale_color_manual(values=(brewer.pal(6, "Spectral")))

plot_cell_trajectory(mono_obj, color_by = "culture")+
  scale_color_manual(values=(brewer.pal(7, "Spectral")))


plot_cell_trajectory(g90, color_by = "gestation")+
  scale_color_manual(values=(brewer.pal(6, "Spectral")))

