setwd("/Users/meyennf/Documents/Cambridge/Datasets/huGroundState_2014/Gafni_Anal/")
filenameRAW <- file.choose() 
filenameRAW <- "/Users/meyennf/Documents/Cambridge/Datasets/huGroundState_2014/Gafni_Anal/CGI_methyl.txt"
dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
SampleNum <- length(dataRAW) - 12
data <- dataRAW[,c(2,13:length(dataRAW))]
names(data) <- c("Chromosome","ICM (Guo)","ICM (Smith)","naive H9 (Takashima)",
                   "primed H9 (Takashima)","primed (Gafni)","naive (Gafni) d12", "naive (Gafni) d17")
data <- data[,c(2,3,5,4,6,7,8)]
data <- data[,c(2:length(dataRAW))]
#data <- data[1:2000,]
# Remove Zero Rows
data <- data[apply(data,1,sum) !=0 ,]
# Reduce Row Number
#data <- data[round(runif(length(data[,1]),0,.52)),]

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
  ggtitle("prcomp")
  

library(FactoMineR)
pca2 = PCA(t(data), graph = T)
pca2$eig
pca2$var$coord
head(pca2$ind$coord)

#VISUALIZE
library(ggplot2)
scores2 = as.data.frame(pca2$ind$coord)

qplot(data = scores2,x = Dim.1, y = Dim.2,xlim=c(-150,150), ylim=c(-150,150),label = rownames(scores1)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "red", alpha = 0.8, size = 4) +
  theme_bw() +
  ggtitle("PCA")

