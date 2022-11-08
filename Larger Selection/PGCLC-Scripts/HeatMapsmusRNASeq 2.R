

library(gplots)
library(RColorBrewer)

#MOUSE
setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/20160621_DiffExp/")
fileName <- "Exp_mEpiLC d2 vs PGCLC d4.txt"
#fileName <- "Exp_mESC naive vs PGCLC d4.txt"


dataRAW <- read.delim(fileName,header=TRUE,skip=0,fill=TRUE,dec = ".")   # Reading in data
data <- dataRAW[,c(1,13:length(dataRAW))]


#MOUSE
DataGroup <- c(rep("naive mESC",3),rep("mEpiLC",3),rep("mPGCLC d4",3),
               rep("mPGCLC d6",3),rep("mPGC E9.5",2),rep("mPGC E11.5",2),
               rep("mPGC E13.5",2),rep("primed mESC",2))

sampleNames <- cbind(DataGroup=DataGroup,
                     Sample=make.unique(c(levels(as.factor(DataGroup)),DataGroup),sep=" ")[-(1:length(levels(as.factor(DataGroup))))])

#sampleNames <- cbind(DataGroup=c(rep("naive mESC",3),rep("mEpiLC",3),rep("mPGCLC",6),
#                                rep("mPGC",6),rep("primed mESC",2)),
#                    Sample=names(data)[-1])


names(data) <- c("Probe",sampleNames[,2])        # assign col (Sample) names

data <- data[,-c(20,21)] # Remove primed ESC
sampleNames <- sampleNames[-c(19,20),] # Remove primed ESC



#data <- na.omit(data)  # Remove Rows with NAs


#data <- data[apply(data[,2:(ncol(data))],1,min)>=0,] #Exclude values samller than 0



mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix

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
my_Col_palette <- colorRampPalette(c("blue", "white","red"))(n = 50)

ColorSetC <- rainbow(length(unique(sampleNames[,1])))
names(ColorSetC) <- unique(sampleNames[,1])
ColColorLabels=ColorSetC[sampleNames[,1]]

ColorSetR <- topo.colors(2*length(unique(groups)))
names(ColorSetR) <- unique(groups)
RowColorLabels=ColorSetR[groups]


# Plotting the heat map
#pdf(paste(fileName,"pdf",sep="."),width=20,height=20,useDingbats=FALSE)
par(cex.lab=1, cex.axis=0.7, cex.sub=0.5)
heatmap.2(mat_data,
          main = paste("Correlation in",fileName), # heat map title
          density="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(8,8),     # widens margins around plot
          col = colorRampPalette(c("yellow"," blue"))(n = 100),       # use on color palette defined earlier
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
#dev.off()


#########################################################
# Export of Data 
#########################################################
output<-cbind(Genes = rownames(mat_data), Cluster = groups)
write.table(output,paste(fileName,"_ClusterGroups.txt",sep=""),sep="\t",row.names = F,col.names = T,quote=T)
remove(output)
#########################################################