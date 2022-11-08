#########################################################
# PCA 
#########################################################

# load libraries
library(ggplot2)
library(grid)

#Copy Data to new DataFrame
data <- dataset

#Remove GeneNames
dataNN <- data[,c(2:ncol(data))]   # dataNN = data with NoNames

# Remove low expressed genes - cutoff are the 10% lowest expressed genes (mean)
dataNN <-  dataNN[(apply(dataNN,1,mean) > quantile(apply(dataNN,1,mean),0.25)),]
hist(apply(dataNN,1,mean))
# Remove Rows with NAs
dataNN <- na.omit(dataNN)
# Remove Zero Rows
# dataNN <- dataNN[apply(dataNN,1,sum) !=0 ,]

#########################################################
# PCA with FACTOMINER
#########################################################
library(FactoMineR)

pcaF = PCA(t(dataNN), graph = F)
  #pcaF$eig
  #head(pcaF$ind$coord)
  #plot(pcaF)

# extract some parts for plotting
PCs <- data.frame(cbind(pcaF$ind$coord[,1],pcaF$ind$coord[,2],pcaF$ind$coord[,3],pcaF$ind$coord[,4]))
names(PCs) <-c("PC1","PC2","PC3","PC4")
rownames(PCs) <- rownames(pcaF$ind$coord)

PCs$Label=sampleNames[,1]
PCs$GroupLabel=sampleNames[,2]

  #ColorSet <- rainbow(length(unique(sampleNames$DataGroup)))
  #names(ColorSet) <- unique(sampleNames$DataGroup)
  #PCs$Color=ColorSet[sampleNames$DataGroup]

qplot(data = PCs,x = PC1, y = PC2,label = Label,cex=1,
      xlab=paste("PC1 (",round(pcaF$eig$`percentage of variance`[1],2),"%)",sep=""),
      ylab=paste("PC2 (",round(pcaF$eig$`percentage of variance`[2],2),"%)",sep=""),
      xlim=c(-150,100), ylim=c(-100,100),
      colour=GroupLabel) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "black", alpha = 0.8, size = 3) +
  theme_bw() +
  ggtitle(paste("PCA",DataSetName))





#########################################################
# PCA with princomp
#########################################################

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
scores = as.data.frame(pca1$x)

qplot(data = scores1,x = PC1, y = PC2,label = rownames(scores1),xlim=c(-10,10), ylim=c(-10,10)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "red", alpha = 0.8, size = 4) +
  theme_bw() +
  ggtitle(paste("PCA","filenameRAW"))



library(FactoMineR)
pca2 = PCA(t(data), graph = F)
plot(pca2)
pca2$eig
pca2$var$coord
head(pca2$ind$coord,n=40)


#VISUALIZE
library(ggplot2)
scores2 = as.data.frame(pca2$ind$coord)

qplot(data = scores2,x = Dim.1, y = Dim.2,label = rownames(scores2)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "black", alpha = 0.8, size = 3) +
  theme_bw() +
  ggtitle(paste("PCA","filenameRAW"))
