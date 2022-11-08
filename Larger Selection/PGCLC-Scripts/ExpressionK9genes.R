library("gplots")


setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/Figures/Figure3/DNMT_Expression/")

#Human DataSets 
fileName <- "K9Genes_HumanExp.csv"

#Mouse DataSets
fileName <- "K9Genes_MouseExp.csv"

mRNAValues <- read.csv(fileName,sep=",",header=TRUE,skip=0,fill=TRUE,dec = ".")


rownames(mRNAValues) <-   mRNAValues$Sample
mRNAValues <- mRNAValues[-1]


data <- t(as.matrix(mRNAValues))
heatmap.2(data,trace="none", density="none",
          main=fileName,
          scale="row",
          Colv = F,Rowv = F,
          dendrogram="none",
          srtCol=45,cexCol=0.7,offsetCol=-0.5,
          col = colorRampPalette(c("  blue","red"))(25))





