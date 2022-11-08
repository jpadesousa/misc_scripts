
setwd("~/Documents/Cambridge/DataSeq/201510_Islets/RNASeq_Islets")
filenameRAW <- "RNASeq_Islets_OBExpValues.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
dataRAW <- dataRAW[,c(1,13,14,15)]
names(dataRAW) <- c("Gene","OBOB","OBOB Keto","wildtype")

SampleNum <- 3

data <- dataRAW[-1]

# Remove Zero Rows
data <- data[apply(data,1,sum) !=0 ,]

# Pricipal Components Analysis
fit <- princomp(data, cor=F)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
#fit$scores # the principal components

pca1 = prcomp(t(data), scale. = T)
pca1$sdev
head(pca1$rotation)
head(pca1$x)

#VISUALIZE
library(ggplot2)
scores1 = as.data.frame(pca1$x)
,xlim=c(-150,150), ylim=c(-150,150)
qplot(data = scores1,x = PC1, y = PC2,label = rownames(scores1)) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "red", alpha = 0.8, size = 4) +
  theme_bw() +
  ggtitle("PCA")
  

library(FactoMineR)
pca2 = PCA(t(data), graph = F)
pca2$eig
#pca2$var$coord
head(pca2$ind$coord)

#VISUALIZE
library(ggplot2)
scores2 = as.data.frame(pca2$ind$coord)

pdf(paste("RNASeq_PCA.pdf"),width=5,height=5,useDingbats=FALSE)
qplot(data = scores2,x = Dim.1, y = Dim.2,label = rownames(scores1),xlim=c(-250,150), ylim=c(-150,150),
      xlab = paste("PC1:",round(pca2$eig$`percentage of variance`[1],2),"%"), 
      ylab = paste("PC2:",round(pca2$eig$`percentage of variance`[2],2),"%")
      ) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  geom_text(colour = "red", alpha = 1, size = 4) +
  theme_bw() +
  ggtitle("Gene expression")
dev.off()
