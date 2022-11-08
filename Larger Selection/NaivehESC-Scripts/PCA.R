require(grid)

setwd("~/Desktop/RNASeq/")


filenameRAW <- "hPGCLC_RNASeqQuantAll.txt"
filenameRAW <- "hPGCLC_RNASeqQuantMatchedDist.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

data <- dataRAW[,c(1,13:length(dataRAW))]

# SamplesNames <- read.table("mPGCLC_Samples.txt",as.is=T)[,1]

# names(dataRAW) <- SamplesNames

# data <- dataRAW


# Remove Rows with NAs
data <- na.omit(data)
# Remove Zero Rows
data <- data[apply(data,1,sum) !=0 ,]


# Pricipal Components Analysis
fit <- princomp(data, cor=F)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components

pca1 = prcomp(t(data), scale. = T)
pca1$sdev
head(pca1$rotation)
head(pca1$x)

#VISUALIZE
library(ggplot2)
scores1 = as.data.frame(pca1$x)

qplot(data = scores1,x = PC1, y = PC2,label = rownames(scores1),xlim=c(-150,150), ylim=c(-150,150)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "red", alpha = 0.8, size = 4) +
  theme_bw() +
  ggtitle(paste("PCA",filenameRAW))



library(FactoMineR)
pca2 = PCA(t(data), graph = F)
plot(pca2)
pca2$eig
pca2$var$coord
head(pca2$ind$coord,n=40)


#VISUALIZE
library(ggplot2)
scores2 = as.data.frame(pca2$ind$coord)

qplot(data = scores2,x = Dim.1, y = Dim.2,xlim=c(-250,250), ylim=c(-250,250),label = rownames(scores1)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "red", alpha = 0.8, size = 4) +
  theme_bw() +
  ggtitle(paste("PCA",filenameRAW))
