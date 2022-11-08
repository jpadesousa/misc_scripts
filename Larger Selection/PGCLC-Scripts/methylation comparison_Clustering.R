#########################################################
# K-means clustering of repeat-free tiles.
#########################################################



#########################################################

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

#########################################################


setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/humBSSeq/methylation comparison/")

#Human DataSets 
fileName <- "NoRepeatsHUM_MoreMore.txt"

#Keep a record of random clustering starting variable
randomVariable <- .Random.seed
MethValues <- read.delim(fileName,header=TRUE,skip=0,fill=TRUE,dec = ".")   # Reading in data

MethValues <- MethValues[,c(1,13:ncol(MethValues))]



# Prepare Data
mat_data <- data.matrix((MethValues[,2:ncol(MethValues)][-9]))  # transform into a matrix
mat_data <- na.omit(mat_data) # listwise deletion of missing
#mat_data <- t(mat_data) # transpose matrix
#mat_data <- scale(mat_data) # standardize variables; methylation is between 1 and 100
head(mat_data)



# Determine number of clusters
wssplot(mat_data)

# 7 cluster solution
numCluster <- 7


# K-Means Cluster Analysis


#set.seed(123)
randomVariable <- cbind(randomVariable,.Random.seed)
randomVariableTime <- Sys.time()
colnames(randomVariable)[ncol(randomVariable)] <- randomVariableTime

#.Random.seed <- attr(fit, "seed")
fit <- kmeans(mat_data, numCluster) 
attr(fit,"seed") <- x

# get cluster means 
aggregate(mat_data,by=list(fit$cluster),FUN=mean)

# append cluster assignment
plotdata <- data.matrix(data.frame(mat_data, fit$cluster))

# Use only part of data
plotdata <- plotdata[1:500,]
head(plotdata,20)

# OrderData by Cluster
plotdata <- plotdata[order(plotdata[,ncol(plotdata)]),]
head(plotdata,20)

rownames(plotdata) <- rep("",nrow(plotdata))

library(gplots)

sampleNames <- cbind(DataGroup=c(rep("hESC",2),rep("hEpiLC",4),rep("hPGCLC",2),
                                 rep("hPGC",3)),
                     Sample=colnames(plotdata)[-12])

ColorSetC <- rainbow(length(unique(sampleNames[,1])))
names(ColorSetC) <- unique(sampleNames[,1])
ColColorLabels=ColorSetC[sampleNames[,1]]

pdf(paste(fileName,"pdf",sep="."),width=20,height=20,useDingbats=FALSE)
par(cex.lab=1, cex.axis=0.7, cex.sub=0.5)

heatmap.2(plotdata[,1:(ncol(plotdata)-1)],
          main = "K-means clustering of repeat-free tiles", # heat map title
          density="none",            # turns off density plot inside color legend
          trace="none",              # turns off trace lines inside the heat map
          margins = c(8,8),         # widens margins around plot
          col = colorRampPalette(c(" dark blue","light yellow", " red"))(n = 100),  # use an color palette defined earlier
          dendrogram = "none",     # no dendrogram
          Rowv=F,Colv=F,
          RowSideColors=rainbow(n = numCluster)[plotdata[,ncol(plotdata)]],
          ColSideColors = ColColorLabels,
          scale = "none",
          srtCol=45,cexCol=0.7,offsetCol=-0.5)            


legend("top", legend = unique(sampleNames[,1]), 
       col = ColorSetC, lty= 1, lwd = 5, cex=0.7)          
dev.off()
