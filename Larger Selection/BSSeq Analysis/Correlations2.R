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
          main=SampleName,
          scale="none",
          labRow="",
          dendrogram="column",
          hclustfun=function(x) hclust(x,method="complete"),
          distfun=function(x) dist(x,method ="euclidean"),
          col = colorRampPalette(c("black", "red"))(100),srtCol=0,
          key.xlab="%CpG methylation", key.title=""
)