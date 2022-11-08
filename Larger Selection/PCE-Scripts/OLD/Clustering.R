#########################################################
# Hierachical Clustering
#########################################################

# Irie: Before clustering and principal component analysis, the transcripts with the 10% lowest average expression were removed, 
# and the gene expression data matrix was centered and scaled. Time points were clustered hierarchically using Ward’s method 
# and the Pearson correlation coefficient (1-c) as a distance measure. Principal component analysis was performed by a singular 
# value decomposition (SVD) of the center-scaled gene expression data matrix.

# Sasaki: For the analysis of differentially expressed genes (DEGs), UHC analysis (hclust function with Euclidian distances and 
# Ward distance functions in R3.1.1), and PCA (prcomp function in R3.1.1)

#########################################################
# Reading in data and transform it into matrix format
#########################################################

#Copy Data to new DataFrame
data <- datasetMus         #Copy Data to new DataFrame
Samples <- SamplesMus      #Copy SampleNames

data <- datasetHum         #Copy Data to new DataFrame
Samples <- SamplesHum      #Copy SampleNames


# rownames(data) <- make.unique(data[,1])            # assign row names
# data <- data [-1]
# colnames(data) <- Samples[,1]
# data <-  data[(apply(data,1,mean) > quantile(apply(data,1,mean),0.10)),]
# 
#  


mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- data[,1]                # assign row names
 
#colnames(mat_data) <- sampleNames[,2]
#colnames(mat_data) <- sampleNames[,1]

# Remove Rows with NAs
mat_data <- na.omit(mat_data)

# Remove low expressed genes
mat_data <-  mat_data[(apply(mat_data,1,mean) > quantile(apply(mat_data,1,mean),0.10)),]

# Standardize variables
mat_data <- scale(mat_data) 

#Random selection of probes to speed process up
#sampling <- 1000
#if (nrow(dataNN) < sampling) { sampling <- nrow(dataNN)}
#mat_data <- mat_data[sample(nrow(mat_data), sampling), ]

#########################################################
# Clustering with pvclust
#########################################################

#install.packages("pvclust")
library("pvclust")

distfunc <- function(x) as.dist((1-cor(x))/2)  # Pearson
result <- pvclust(mat_data, method.dist=distfunc, method.hclust="complete", nboot=1000,parallel=T)

#result <- pvclust(mat_data, method.dist="euclidean", method.hclust="complete", nboot=1000,parallel=T)

plot(result)
pvrect(result, alpha=0.95)

#########################################################
# Clustering with hclust - Euclidian distances and Ward distance func
#########################################################

#distfunc <- function(x) dist(t(x),method="euclidean")
distfunc <- function(x) as.dist((1-cor(x))/2)  # Pearson

result <- hclust(distfunc(mat_data),method="ward.D")
#result <- hclust(distfunc(mat_data),method="complete")

plot(result)


#########################################################
# K-means clustering
#########################################################
inv_data <- t(mat_data)
# Determine number of clusters
wss <- (nrow(inv_data)-1)*sum(apply(inv_data,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(inv_data, 
                                     centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")



#########################################################
# SC3 clustering
#########################################################
#biocLite("SC3")
library(SC3)

sc3(data,3:7,log.scale=F,show.original.labels=T,n.cores=8)

sc3(data, ks = 3:7, cell.filter = FALSE, cell.filter.genes = 2000,
    gene.filter = TRUE, gene.filter.fraction = 0.06, log.scale = TRUE,
    d.region.min = 0.04, d.region.max = 0.07, interactivity = TRUE,
    show.original.labels = FALSE, svm = FALSE, svm.num.cells = NA,
    n.cores = NA, seed = 1)

