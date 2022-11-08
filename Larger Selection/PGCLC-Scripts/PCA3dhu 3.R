#########################################################
# Import of RNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/humRNASeq/20160624RPKM/")

fileName <- "Expression_humPGCLC_RNASeq.txt"
dataset <- read.csv(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

dataset <- dataset[,c(1,13:ncol(dataset))]

DataGroup <- c(rep("naive hESC",3),rep("primed hESC",3),rep("hPGCLC",5),rep("hPGC",6),rep("hSoma",3))

#sampleNames <- make.unique(c(levels(as.factor(DataGroup)),DataGroup),sep=" ")[-(1:length

sampleNames <- c("naive hESC 1","naive hESC 2","naive hESC 3",
                 "primed hESC 1","primed hESC 2","primed hESC 3",
                 "hPGCLC d4 1","hPGCLC d4 2","hPGCLC d5",
                 "hPGCLC d8","hPGCLC d12","hPGC Wk5.5",
                 "hPGC Wk7F 1","hPGC Wk7F 2","hPGC Wk7F 3",
                 "hPGC Wk9F 1","hPGC Wk9F 2",
                 "hSoma 1","hSoma 2","hSoma 3")


#Samples <- cbind(DataGroup,names(dataset)[-1])
Samples <- cbind(DataGroup,sampleNames)
names(dataset) <- c("Probe",sampleNames)
data <- dataset
#########################################################

data <- datasetMus         #Copy Data to new DataFrame
Samples <- SamplesMus      #Copy SampleNames

#data <- datasetHum         #Copy Data to new DataFrame
#Samples <- SamplesHum      #Copy SampleNames

dataNN <- data[,c(2:ncol(data))]   # Remove GeneNames: dataNN = data with NoNames #
names(dataNN)
#Remove or Select DataSets
  SELECT <- c(13:15,16:18, 19:22 )
  dataNN <- dataNN[,-SELECT]
  Samples <- Samples[-SELECT,]


#hist(apply(dataNN,1,mean)) #Distribution of expression values

dataNN <-  dataNN[(apply(dataNN,1,mean) > quantile(apply(dataNN,1,mean),0.10)),] # Remove low expressed genes

#########################################################
# PCA with FACTOMINER
#########################################################
library(FactoMineR)

dataNN <- na.omit(dataNN)       # Remove Rows with NAs

pcaF = PCA(t(dataNN), graph = F)
#pcaF$eig
#head(pcaF$ind$coord)


pcaF$ind$coord[,1] <- -pcaF$ind$coord[,1]
#pcaF$ind$coord[,2] <- -pcaF$ind$coord[,2]
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
                       main="PCA human RNASeq",
                       xlab="PC1",
                       ylab="PC2",
                       zlab="PC3")
  s3d.coords <- s3d$xyz.convert(PC1, PC2, PC3) # convert 3D coords to 2D projection
  text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
       labels=row.names(PCs),               # text to plot
       cex=.5, pos=4)           # shrink text 50% and place to right of points)
})

