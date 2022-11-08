setwd("~/Documents/Cambridge/Collaborations/Austin Smith/2015Oct/Data")
filenameRAW <- file.choose() 
SampleName <- strsplit(strsplit(filenameRAW, "/")[[1]][10],".txt")[[1]]

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".") #,nrows=50
data <- dataRAW[,c(2,13:length(dataRAW))]
SampleNum <- length(data)-1

#data <- data[1:100,]

#Exclude non-quantified samples
data <- data[apply(data[,2:(1+SampleNum)],1,min)>=0,]


#data$Range <- abs(apply(data[,c(3,4)],1,mean) - apply(data[,c(5,6)],1,mean))
#data$Range <- apply(data[,5:6],1,max) - apply(data[,2:3],1,min)

#Random selection of probes
sampling <- 10000
if (nrow(data)< 10000) { sampling <- nrow(data)}
data2Use <- data[sample(nrow(data), sampling), ]

#data2Use <- data[data$Range>80,]
#data2Use <- data[round(runif(length(data[,1]),0,.52)),]
#data_matrix <- data.matrix(data[,2:(1+SampleNum)])

data_matrix <- data.matrix(data2Use[,2:(1+SampleNum)])

#head(data_matrix)

library("gplots")
heatmap.2(data_matrix,trace="none", density="none",
          scale="row", #note this is important because it allows you to use RowZ-score, or actual values
          labRow="",
          dendrogram="col",
          hclustfun=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          col = colorRampPalette(c("black", "red"))(100),srtCol=0
)
heatmap.2(data_matrix,trace="none", density="none",
          main=SampleName,
          scale="none",
          labRow="",
          dendrogram="column",
          hclustfun=function(x) hclust(x,method="complete"),
          distfun=function(x) dist(x,method ="euclidean"),
          col = colorRampPalette(c("black", "red"))(100),srtCol=0,
          key.xlab="%CpG methylation", key.title=""
)
library("pvclust")

pvclust
result <- pvclust(data_matrix, method.dist="euclidean", method.hclust="complete", nboot=1000)
plot(result)
pvrect(result, alpha=0.95)



mydata <- na.omit(data_matrix) # listwise deletion of missing
mydata <- scale(mydata) # standardize variables

# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(mydata, 2) # 2 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)

# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(mydata, method.hclust="ward.D",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 2)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit$cluster)

