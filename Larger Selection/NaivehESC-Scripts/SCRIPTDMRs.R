library(gplots)
library(RColorBrewer)

setwd("/Users/meyennf/Documents/Cambridge/Projects/hESC/Guo_2/20170505/")



#DMRs
filenameRAW <- "DMR_meth.txt"

Data <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
Data <- Data[,c(1,13:ncol(Data))]

rownames(Data) <- Data[,1]
Data <- t(Data[,-1])
rownames(Data) <- DataGroups <- c(rep("Sperm",3),rep("Oocyte",1),rep("ICM",2),
                                  rep("Shef6-EOS",3),rep("cR-Shef6-EOS_p5",3),rep("cR-Shef6-EOS_p9",3),rep("cR-Shef6-EOS_LT",2),
                                  rep("H9-NK2",3),rep("H9-NK2-reset",3),rep("HNES1",3),rep("HNES1-primed",3),
                                  rep("cR-H9.1",3),rep("cR-H9.2",3),rep("cR-S6.2",3))


#DMRs
AggData <- aggregate(Data, by = list(DataGroups),
          mean, na.rm = TRUE) 

AggData <- AggData[c(14,12,11,13,5,6,4,7,8,9,10,1,2,3),]
rownames(AggData) <- AggData[,1]

AggData <- t(AggData[,-1])
mat_data <-data.matrix(AggData)

#rowMeans(Data)

#mat_data <- data.matrix(AggData[,2:ncol(AggData)])  # transform column 2-5 into a matrix
#rownames(mat_data) <- AggData[,1]                # assign row (Gene) names

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




