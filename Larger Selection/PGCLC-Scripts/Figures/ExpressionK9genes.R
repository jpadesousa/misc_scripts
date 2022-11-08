library("gplots")
setwd("/Users/meyennf/Desktop//Figure4/K9Expression/")


#Human DataSets 
fileName <- "K9Genes_HumanExp.csv"

#Mouse DataSets
fileName <- "K9Genes_MouseExp.csv"

mRNAValues <- read.csv(fileName,sep=",",header=TRUE,skip=0,fill=TRUE,dec = ".")

rownames(mRNAValues) <-   mRNAValues$Probe
mRNAValues <- mRNAValues[-1]

#pdf (paste(fileName,"K9.pdf"), width = 10, height = 10, useDingbats=FALSE,colormodel='rgb')

par(cex.lab=1, cex.axis=0.7, cex.sub=0.5, mar=c(5,5,3,3))


data <- t(as.matrix(mRNAValues))
heatmap.2(data,trace="none", density="none",
          main=fileName,
          scale="none",
          Colv = F,Rowv = F,
          dendrogram="none",
          srtCol=45,cexCol=0.7,offsetCol=-0.5,
          col = colorRampPalette(c("white","dark blue"))(n = 25)
)

#dev.off()

