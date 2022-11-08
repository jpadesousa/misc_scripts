#########################################################
# Import of RNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/")

DataSetName <- "musPGCLC_RNASeq_OUT_20160614.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

dataset <- datasetInput[,c(1,13:ncol(datasetInput))]
remove(datasetInput)

dataNames <- cbind(c(rep("Seisenberger",2),rep("Yamaguchi",4),rep("ESC",3),rep("EpiLC",3),rep("PGCLC",8)),names(dataset)[-1])

#########################################################

#Copy Data to new DataFrame
data <- dataset
sampleNames <- dataNames

#Remove GeneNames
dataNN <- data[,c(2:ncol(data))]   # dataNN = data with NoNames


#Remove Merged DataSets
#data <- data[-(2:5)]
#sampleNames <- dataNames[-(1:4),]

#Reduce DataSet Number
rem <- c(3:6,9,12,15,16,19,20)
dataNN <- dataNN[-rem]
sampleNames <- sampleNames[-rem,]


#Distribution of expression values
hist(apply(dataNN,1,mean))

# Remove low expressed genes - cutoff are the 10% lowest expressed genes (mean)
dataNN <-  dataNN[(apply(dataNN,1,mean) > quantile(apply(dataNN,1,mean),0.1)),]
hist(apply(dataNN,1,mean))

# Remove Rows with NAs
dataNN <- na.omit(dataNN)

# Remove Zero Rows
dataNN <- dataNN[apply(dataNN,1,sum) !=0 ,]

#########################################################
# PCA with FACTOMINER
#########################################################
library(FactoMineR)

pcaF = PCA(t(dataNN), graph = F)
#pcaF$eig
#head(pcaF$ind$coord)
plot(pcaF)

# extract some parts for plotting
PCs <- data.frame(cbind(pcaF$ind$coord[,1],pcaF$ind$coord[,2],pcaF$ind$coord[,3],pcaF$ind$coord[,4]))
names(PCs) <-c("PC1","PC2","PC3","PC4")
rownames(PCs) <- rownames(pcaF$ind$coord)

# 3D Scatterplot
library(scatterplot3d)
with(PCs, {
  s3d <- scatterplot3d(PC1, PC2, PC3,        # x y and z axis
                       highlight.3d=TRUE,pch=18,  cex.symbols=2,   # filled  circles
                       type="h",                    # vertical lines to the x-y plane
                       main="PCA mouse RNASeq",
                       xlab="PC1",
                       ylab="PC2",
                       zlab="PC3")
  s3d.coords <- s3d$xyz.convert(PC1, PC2, PC3) # convert 3D coords to 2D projection
  text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
       labels=row.names(PCs),               # text to plot
       cex=.5, pos=4)           # shrink text 50% and place to right of points)
})


