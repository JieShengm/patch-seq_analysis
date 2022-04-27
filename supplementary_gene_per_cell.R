setwd("~/Desktop/Xinyu_Patchseq/")
ref = read.table("./reference_genepercellcounts.txt")
human = read.table("./human_genepercellcounts.txt",row.names = 1,skip =  1)
monkey = read.table("./monkey_genepercellcounts.txt",row.names = 1,skip =  1)
mouse = read.table("./reference_data/mouse/mouse_genepercellcounts.txt", row.names = 1, skip = 1)
df = data.frame(count = c(ref$x, human$V2, monkey$V2, mouse$V2), 
                group = c(rep("BICCN",dim(ref)[1]),
                          rep("human",dim(human)[1]),
                          rep("monkey",dim(monkey)[1]),
                          rep("mouse", dim(mouse)[1])))
library(ggplot2)
library(scales)

pdf("nGene_per_cell_w_mouse.pdf")
ggplot(df,aes(x=count,color = group,fill=group))+
  geom_histogram(aes(y=1000*..density..),
                 alpha=0.4,position='identity',binwidth=1000)+
  scale_y_continuous(labels = percent_format(), name = "percent")+
  xlab("The number of genes per cell")+
  scale_color_manual(values=c("#999999", "#A93226", "#2980B9","#F5B041"))+
  scale_fill_manual(values=c("#999999", "#A93226", "#2980B9","#F5B041"))+ 
  theme_classic()+
  scale_x_continuous(limits = c(0,24000))
dev.off()  

hist(monkey$V2,breaks = seq(0,24000,1000))
h = hist(monkey$V2,breaks = seq(0,24000,1000)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)

hist(human$V2,breaks = seq(0,24000,1000))
h = hist(human$V2,breaks = seq(0,24000,1000)) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)
