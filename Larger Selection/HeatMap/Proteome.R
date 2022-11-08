library("gplots")

setwd("~/Dropbox/Manuscripts/Serum2i/Figures/Figure 4/Proteome_Uhrf1/")
fileName <- "Proteome.txt"
Proteome <- read.table(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

#Proteome <- cbind(Proteome[1],rep(0,nrow(Proteome)),Proteome[,2:length(Proteome)])

names(Proteome) <-c("Gene",4,8,16,24,32,72)

Proteome[1]

Samples <- c(1,2,3,4,17,18,21)

data <- as.matrix(Proteome[Samples,2:7])
heatmap.2(data,trace="none", density="none",
          scale="none",
          Colv = F,Rowv = F,
          dendrogram="none",
          labRow=Proteome[Samples,1],
          
          col = colorRampPalette(c(" blue"," red"))(25),srtCol=0
          )



