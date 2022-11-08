

library(gplots)
library(RColorBrewer)

#HUMAN
setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/humRNASeq/20160621_DiffExp/")
fileName <- "Exp_hESC primed hPGCLC.txt"
#fileName <- "Exp_hESC primed naive.txt"




dataRAW <- read.delim(fileName,header=TRUE,skip=0,fill=TRUE,dec = ".")   # Reading in data
data <- dataRAW[,c(1,13:length(dataRAW))]

#HUMAN
sampleNames <- cbind(DataGroup=c(rep("naive hESC",3),
                                 rep("primed hESC",3),
                                 rep("hPGCLC",5),
                                 rep("hPGC",6),
                                 rep("hSoma",3)),
                     Sample=c("naive 1","naive 2","naive 3",
                              "primed 1","primed 2","primed 3",
                              "hPGCLC d4","hPGCLC d4 2","hPGCLC d5","hPGCLC d8","hPGCLC d12",
                              "hPGC Wk5.5","hPGC Wk7F 1","hPGC Wk7F 2","hPGC Wk7F 3","hPGC Wk9F 1","hPGC Wk9F 2",
                              "hSoma 1","hSoma 2","hSoma 3"))

names(data) <- c("Probe",sampleNames[,2])        # assign col (Sample) names

data <- na.omit(data) # Remove Rows with NAs

data[,8] <- apply (data[,c(8,9)],1,mean) # Average of hPGCLC d4
data <- data[,-c(9,17:21)] # Remove 2nd hPGCLC d4, hSoma, hPGC Wk9
sampleNames <- sampleNames[-c(8,16:20),] # Remove 2nd hPGCLC d4, hSoma, hPGC Wk9

#Exclude non-quantified samples
#data <- data[apply(data[,2:(ncol(data))],1,min)>=0,]

mat_data <- data.matrix(data[,2:ncol(data)])  # transform into a matrix

rownames(mat_data) <- data[,1]                # assign row (Gene) names


# Clustering

#hclustfunc <- function(x) hclust(x, method="complete")
hclustfunc <- function(x) hclust(x, method="ward.D")

distfunc <- function(x) dist(x,method="euclidean")
#distfunc <- function(x) as.dist((1-cor(t(x)))/2) # Pearson

clusteringResult <- hclustfunc(distfunc(mat_data))
#result <- hclust(distfunc(mat_data),method="complete")

clusterNum <- 4
groups<-cutree(clusteringResult, k=clusterNum)

# Customizing and plotting the heat map

# creates a own color palette from red to green
my_Col_palette <- colorRampPalette(c("blue", "white","red"))(n = 299)

ColorSetC <- rainbow(length(unique(sampleNames[,1])))
names(ColorSetC) <- unique(sampleNames[,1])
ColColorLabels=ColorSetC[sampleNames[,1]]

ColorSetR <- topo.colors(length(unique(groups)))
names(ColorSetR) <- unique(groups)
RowColorLabels=ColorSetR[groups]


# Plotting the heat map
pdf(paste(fileName,"pdf",sep="."),width=20,height=20,useDingbats=FALSE)
par(cex.lab=1, cex.axis=0.7, cex.sub=0.5)
heatmap.2(mat_data,
          main = paste("Correlation in",fileName), # heat map title
          density="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(8,8),     # widens margins around plot
          col = colorRampPalette(c("white","dark blue"))(n = 100),       # use on color palette defined earlier
          dendrogram = "row",     # only draw a row or col dendrogram
          scale = "row",
          Colv=F,Rowv=T,
          RowSideColors = RowColorLabels,
          ColSideColors = ColColorLabels,
          hclustfun = hclustfunc,distfun = distfunc,
          srtCol=45,cexCol=0.7,offsetCol=-0.5)            

legend("top", legend = unique(sampleNames[,1]), 
       col = ColorSetC, lty= 1, lwd = 5, cex=0.7)   
legend("bottomleft", legend = unique(groups), 
       col = ColorSetR, lty= 1, lwd = 5, cex=0.7)          
dev.off()




#########################################################
# Export of Data 
#########################################################
output<-cbind(Genes = rownames(mat_data), Cluster = groups)
write.table(output,paste(fileName,"_ClusterGroups.txt",sep=""),sep="\t",row.names = F,col.names = T,quote=T)
remove(output)
#########################################################