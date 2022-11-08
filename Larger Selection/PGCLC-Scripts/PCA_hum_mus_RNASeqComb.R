#########################################################
# Import of huRNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/humRNASeq/20160624_2/")

fileName <- "Expression_humPGCLC_RNASeq.txt"
dataset <- read.csv(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

datasetHum <- dataset[,c(1,13:ncol(dataset))]

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
SamplesHum <- cbind(DataGroup,sampleNames)
names(datasetHum) <- c("Probe",sampleNames)
remove(DataGroup,fileName,sampleNames,dataset)
#########################################################
#Remove Samples
datasetHum <- datasetHum[,-c(5,6,7, 17,18, 19,20,21 )]

#########################################################
# Import of mus RNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/20160624/")

fileName <- "Expression_musPGCLC_RNASeq.txt"
dataset <- read.csv(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

datasetMus <- dataset[,c(1,13:ncol(dataset))]

DataGroup <- c(rep("naive mESC",3),rep("mEpiLC",3),rep("mPGCLC d4",3),
               rep("mPGCLC d6",3),rep("mPGC E9.5",2),rep("mPGC E11.5",2),
               rep("mPGC E13.5",2),rep("primed mESC",2),rep("naive mESC",2))

sampleNames <- make.unique(c(levels(as.factor(DataGroup)),DataGroup),sep=" ")[-(1:length(levels(as.factor(DataGroup))))]

#Samples <- cbind(DataGroup,names(dataset)[-1])

SamplesMus <- cbind(DataGroup,sampleNames)
names(datasetMus) <- c("Probe",sampleNames)
remove(DataGroup,fileName,sampleNames,dataset)
#########################################################

#Remove Samples
datasetMus <- datasetMus[,-c(20,21, 22,23)]
c(13:15,16:18, 19:22 )
#########################################################


match(toupper(datasetMus$Probe),toupper(datasetHum$Probe)) -> pos

datasetComb <- cbind(datasetHum[pos,],datasetMus[-1])
datasetComb <- na.omit(datasetComb) # Remove Rows with NAs
rownames(datasetComb)  <- make.unique(datasetComb[,1])
dataNN <- datasetComb[,c(2:ncol(datasetComb))]   # dataNN = data with NoNames





#Distribution of expression values
hist(apply(dataNN,1,mean))

# Remove low expressed genes - cutoff are the 10% lowest expressed genes (mean)
dataNN <-  dataNN[(apply(dataNN,1,mean) > quantile(apply(dataNN,1,mean),0.25)),]
hist(apply(dataNN,1,mean))


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

cbind (rownames(PCs),pcaF$call$X$TFAP2C)

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

library(ggplot2)
PCs$Label=names(dataNN)  

qplot(data = PCs,x = PC2, y = PC3,label = Label,cex=1,
      xlab=paste("PC1 (",round(pcaF$eig$`percentage of variance`[1],2),"%)",sep=""),
      ylab=paste("PC2 (",round(pcaF$eig$`percentage of variance`[2],2),"%)",sep=""),
      xlim=c(-150,100), ylim=c(-100,100)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "black", alpha = 0.8, size = 3) +
  theme_bw() +
  ggtitle(paste("PCA"))
