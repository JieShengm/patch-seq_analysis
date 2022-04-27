setwd("~/Desktop/Xinyu_Patchseq/")

library(reshape2)
library(ManiNetCluster)
library(plyr)
library(RColorBrewer)
Dim_red = function(edata,gdata,method,cellnames,d,k_NN,k_medoid){
  #e-feature
  X = apply(edata,2,scale)
  #g-feature
  Y=t(log2(gdata+1))
  #Dim_red
  XY_corr=Correspondence(matrix=diag(nrow(X)))
  df=ManiNetCluster(X,Y,nameX='Ephys',nameY='Expr',corr=XY_corr,d=d,
                    method=method,k_NN=k_NN,k_medoids=k_medoid)
  df$cellnames = rep(unlist(cellnames),2)
  return(df[,-1])
}

# load data & gene markers
gdata = read.csv("./human_gdata.csv",row.names = 1)
colnames(gdata) = gsub("X","",colnames(gdata))
edata = read.csv("./human_edata.csv",row.names = 1)

#feature selection elec
edata$Days.in.vitro..DIV.=NULL
edata$Gestation.Day..GD.=NULL
edata$UMAP1=NULL
edata$UMAP2=NULL

#feature selection gene

#gene_type$Cluster = sapply(gene_type$Cluster,substr,start=1,stop=2)
#genenames = unique(gene_type$Gene[gene_type$Cluster %in% c("Ex","In")])
#genenames = na.omit(rownames(gdata)[match(genenames,toupper(rownames(gdata)))]) #1298
#gdata = gdata[intersect(genenames,rownames(gdata)),cellnames]

n = nrow(edata)
cellnames = rownames(edata)

### MA ###
NMA_res = Dim_red(edata,gdata,method = 'nonlinear manifold aln',cellnames =cellnames,d=3L,k_NN=2L,k_medoid=5L)
NMA_res_e = NMA_res[NMA_res$data=="Ephys",]
NMA_res_t= NMA_res[NMA_res$data=="Expr",]

write.csv(NMA_res_e[,2:4],"human_NMA_e.csv")
write.csv(NMA_res_t[,2:4],"human_NMA_t.csv")


#Plot result
library(plot3D)
#Figure S1
#NMA
NMA_plot = -rbind(apply(NMA_res[1:n,2:4],2,scale),apply(NMA_res[(n+1):(2*n),2:4],2,scale))
NMA_plot[1:n,1] = NMA_plot[1:n,1] + 0.5
points3D(x=NMA_plot[,1], y=NMA_plot[,2], z=NMA_plot[,3],pch = 19,cex=0.5,bty="g",ticktype = "detailed", theta = 40, phi = 10,
         xlab = "",ylab = "",zlab = "",box=T, axes=F,
         colvar = as.numeric(mapvalues(NMA_res$data,names(table(NMA_res$data)),1:2)),col = brewer.pal(6,"Spectral")[c(1,6)],
         colkey = F)
