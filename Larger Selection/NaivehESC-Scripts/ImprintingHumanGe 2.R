
library(gplots)
library(RColorBrewer)

setwd("/Users/meyennf/Documents/Cambridge/Projects/hESC/Guo_2/Imprints/")

filenameRAW <- "Export_ImprintedDMR_Expanded.txt"
#filenameRAW <- "Export_ImprintedDMR_Grouped.txt"
#filenameRAW <- "Export_KnownImprintedDMR_Expanded.txt"
#filenameRAW <- "Export_KnownImprintedDMR_Grouped.txt"
  filenameRAW <- "Export_ImprintedDMRComp20170413.txt"


Imprints <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")

Imprints <- Imprints[,c(1,13:ncol(Imprints))]
#names(Imprints) <- c("Imprint","Origin","hSperm","hOocyte","hICM","hESC H9 naive","hEpiLC d1","hEpiLC d2",
#                 "hEpiLC d3","hEpiLC d4","hPGCLC d4","hPGCLC d5","hPGCLC d8","hPGCLC d12","hPGC Wk5.5","hPGC Wk7")

#Imprints <- (  Imprints[,-13])
data <- Imprints
#data <- na.omit(Imprints)
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- data[,1]                # assign row (Gene) names

paste(filenameRAW,"pdf",sep=".")
pdf (paste(filenameRAW,"pdf",sep="."), width = 10, height = 10, useDingbats=FALSE,colormodel='rgb')

par(cex.lab=1, cex.axis=0.7, cex.sub=0.5, mar=c(5,5,3,3))

heatmap.2(mat_data,
          main = paste("Correlation in"), # heat map title
          density="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(8,8),     # widens margins around plot
          col = colorRampPalette(c(" dark blue", " red"))(n = 100),  # use an color palette defined earlier
          dendrogram = "row",     # only draw a row or col dendrogram
 #         scale = "col",
          Colv=F,Rowv=T,
#          RowSideColors = RowColorLabels,
 #         ColSideColors = ColColorLabels,
 #         hclustfun = hclustfunc,distfun = distfunc,
          srtCol=45,cexCol=0.7,offsetCol=-0.5)            





#legend("bottomright", legend = names(Input)[-1],        col = cols, lty= 1, lwd = 5, cex=0.5)          


dev.off()


