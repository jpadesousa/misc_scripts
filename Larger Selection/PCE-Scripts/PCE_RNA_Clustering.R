#####################################################################
# Hierachical Clustering
#####################################################################

#####################################################################
# Reading in data, transform it into matrix format, and preprocess
#####################################################################

#Copy Data to new DataFrame

mat_data <- data.matrix(dataset)                    # transform column 2-5 into a matrix

# Remove Rows with NAs
mat_data <- na.omit(mat_data)

# Remove low expressed genes
#mat_data <-  mat_data[(apply(mat_data,1,mean) > quantile(apply(mat_data,1,mean),0.10)),]

# Remove genes expressed less than 0 log2 RPM in all samples
# mat_data <-  mat_data[(apply(mat_data,1,max) > 0),]

# Remove genes expressed less then 0 log2 RPM in all 4 groups
groupMean <- cbind(apply(mat_data[,c(1:6)],1,mean),apply(mat_data[,c(7:12)],1,mean),
                   apply(mat_data[,c(13:18)],1,mean),apply(mat_data[,c(19:24)],1,mean))
colnames(groupMean) <- Samples[c(1,7,13,18),2]
mat_data <- mat_data[(apply(groupMean,1,max) > 0),]


nrow(dataset)  # Number of Genes analysed
nrow(mat_data) # Number of Genes remaining

# Standardize variables
mat_data <- scale(mat_data) 

#Random selection of probes to speed process up
#sampling <- 1000
#if (nrow(dataNN) < sampling) { sampling <- nrow(dataNN)}
#mat_data <- mat_data[sample(nrow(mat_data), sampling), ]

#########################################################
# Clustering with pvclust
#########################################################

# install.packages("pvclust")
# library("pvclust")

# distfunc <- function(x) as.dist((1-cor(x))/2)  # Pearson
# result <- pvclust(mat_data, method.dist=distfunc, method.hclust="complete", nboot=1000,parallel=T)
# 
# #result <- pvclust(mat_data, method.dist="euclidean", method.hclust="complete", nboot=1000,parallel=T)
# 
# plot(result)
# pvrect(result, alpha=0.95)

#########################################################
# Clustering with hclust - Euclidian distances and Ward distance func
#########################################################

distfunc <- function(x) dist(t(x),method="euclidean")
#distfunc <- function(x) as.dist((1-cor(x))/2)  # Pearson

result <- hclust(distfunc(mat_data),method="ward.D")
#result <- hclust(distfunc(mat_data),method="complete")
#result <- hclust(distfunc(mat_data),method="single")

#########################################################
# Plot and save as PDF
#########################################################
plot(result,hang = -1, cex =1)



#install.packages('dendextend')
library(dendextend)
dend <- as.dendrogram(result)

pdf("Results/RNAseq_Hierachical_Clustering.pdf",width=5,height=5,useDingbats=FALSE)

dend %>% set("leaves_pch", 16) %>% set("leaves_col", Color4Samples[Samples$sampleGroup]) %>% set("labels_cex", 0.7) %>%
  plot(ylab = "Height",main = "Distance:euclidean & Linkage: Ward.D",horiz = FALSE)

dev.off()

remove(dend,groupMean,mat_data,distfunc,result)


