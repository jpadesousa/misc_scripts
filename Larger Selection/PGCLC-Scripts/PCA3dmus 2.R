#########################################################
# Import of RNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/20160624/")

fileName <- "Expression_musPGCLC_RNASeq.txt"
dataset <- read.csv(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

dataset <- dataset[,c(1,13:ncol(dataset))]

DataGroup <- c(rep("naive mESC",3),rep("mEpiLC",3),rep("mPGCLC d4",3),
               rep("mPGCLC d6",3),rep("mPGC E9.5",2),rep("mPGC E11.5",2),
               rep("mPGC E13.5",2),rep("primed mESC",2),rep("naive mESC",2))

sampleNames <- make.unique(c(levels(as.factor(DataGroup)),DataGroup),sep=" ")[-(1:length(levels(as.factor(DataGroup))))]

#Samples <- cbind(DataGroup,names(dataset)[-1])

Samples <- cbind(DataGroup,sampleNames)


#########################################################

#Copy Data to new DataFrame
data <- dataset
names(data) <- c("Gene",sampleNames)


#Remove GeneNames
dataNN <- data[,c(2:ncol(data))]   # dataNN = data with NoNames

#Remove DataSets
# rem <- c(3:6,9,12,15,16,19,20)
# dataNN <- dataNN[-rem]
# Samples <- Samples[-rem,]

#Select DataSets
# SELECT <- c(5,7,12,13,16,17,23,24,25,26,27,29,30)
# dataNN <- dataNN[,SELECT]
# Samples <- Samples[,SELECT]

#Distribution of expression values
hist(apply(dataNN,1,mean))

# Remove low expressed genes - cutoff are the 10% lowest expressed genes (mean)
dataNN <-  dataNN[(apply(dataNN,1,mean) > quantile(apply(dataNN,1,mean),0.25)),]
hist(apply(dataNN,1,mean))

# Remove Rows with NAs
dataNN <- na.omit(dataNN)

# Remove Zero Rows
#dataNN <- dataNN[apply(dataNN,1,sum) !=0 ,]

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


