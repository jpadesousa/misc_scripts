library("gplots")
setwd("/Users/meyennf/Desktop/Exp/")


#fileR <- "hum_SelectedGenesDNMT.txt.csv"

filenames <- list.files(pattern="*.txt.csv", full.names=F)
for (fileR in filenames) {
  mRNAValues <- read.csv(fileR,header=TRUE,skip=0,fill=TRUE,dec = ".")
  rownames(mRNAValues) <-   mRNAValues$Probe
  mRNAValues <- mRNAValues[-1]
  
  data <- t(as.matrix(mRNAValues))
  
  pdf (paste(fileR,".pdf",sep=""), width = 10, height = 10, useDingbats=FALSE,colormodel='rgb')
  par(cex.lab=1, cex.axis=0.7, cex.sub=0.5, mar=c(5,5,3,3))
  heatmap.2(data,trace="none", density="none",
            main=fileR,
            scale="none",
            Colv = F,Rowv = F,
            dendrogram="none",
            srtCol=45,cexCol=0.7,offsetCol=-0.5,
            col = colorRampPalette(c("white"," dark blue"))(n = 25)
  )
  dev.off()
}

