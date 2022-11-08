#####################################################################
# PCA
#####################################################################

#####################################################################
# Reading in data, transform it into matrix format, and preprocess
#####################################################################
source('~/Cambridge/Projects/Collaboration_Wolfrum/2017_PCE/2017October/Scripts/PCE_RNA_PreProcessing.R')


#Copy Data to new DataFrame
#mat_data <- data.matrix(dataset)  # transform into a matrix
mat_data <- dataset  # transform into a matrix

# Remove Rows with NAs
mat_data <- na.omit(mat_data)

# Remove Zero Rows
mat_data <- mat_data[apply(mat_data,1,sum) !=0 ,]

# Remove genes expressed less then 0 log2 RPM in all 4 groups
groupMean <- cbind(apply(mat_data[,c(1:6)],1,mean),apply(mat_data[,c(7:12)],1,mean),
                   apply(mat_data[,c(13:18)],1,mean),apply(mat_data[,c(19:24)],1,mean))
#colnames(groupMean) <- Samples[c(1,7,13,18),2]
mat_data <- mat_data[(apply(groupMean,1,max) > 0),]

nrow(dataset)  # Number of Genes analysed
nrow(mat_data) # Number of Genes remaining

# #####################################################################
# # Pricipal Components Analysis: princomp 
# #####################################################################
# 
# pca.results <- princomp(mat_data,  cor = F, scores = T)
# summary(pca.results)           # print variance accounted for 
# loadings(pca.results)          # pc loadings 
# plot(pca.results,type="lines") # scree plot 
# pca.results$scores             # the principal components


#####################################################################
# Pricipal Components Analysis: prcomp 
#####################################################################
pca.results = prcomp(t(mat_data))

loading <- (eigenvalues <- pca.results$sdev^2) / sum(eigenvalues) * 100
head(eigenvalues)           # Variance data
head(pca.results$rotation)  # Weightings
head(pca.results$x)         # Principal Components 

#VISUALIZE
library(ggplot2)

scores = as.data.frame(pca.results$x)

ggplot(data=scores,aes(x=PC1, y=PC2)) +
  geom_point(pch=16,colour = Color4Samples[Samples$sampleGroup],size = 2) +
  geom_text(label = rownames(scores),colour = "black", alpha = 1, size = 2,hjust=0, vjust=0) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_vline(xintercept = 0, colour = "grey") +
  theme(legend.position="none") +
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste("PC1:",round(loading[1]),"%")) + 
  ylab(paste("PC2:",round(loading[2]),"%")) +
  ggtitle("PCA - prcomp (default settings)")
ggsave("Results/RNASeq_PCA.pdf",width=5,height=5,useDingbats=FALSE)


# #########################################################
# # Pricipal Components Analysis: FactoMineR
# #########################################################  
# # install.packages("FactoMineR")
# library(FactoMineR)
# pca.results = PCA(t(mat_data), graph = F)
# 
# pca.results$eig
# head(pca.results$var$coord)   # Variance data
# head(pca.results$ind$coord)   # Principal Components 
# 
# #VISUALIZE
# library(ggplot2)
# scores = as.data.frame(pca.results$ind$coord)
# 
# qplot(data = scores,x = Dim.1, y = Dim.2,label = rownames(scores),
#       xlab = paste("PC1:",round(pca.results$eig[1,2]),"%"), 
#       ylab = paste("PC2:",round(pca.results$eig[2,2]),"%")) +
#   
#   geom_hline(yintercept = 0, colour = "grey") +
#   geom_vline(xintercept = 0, colour = "grey") +
#   geom_text(colour = "red", alpha = 1, size = 4) +
#   theme_bw() +
#   ggtitle("PCA - princomp")


remove(scores,pca.results,groupMean,mat_data,loading,eigenvalues)




