setwd("~/Desktop/Keystone")
dataRAW <- read.csv(file.choose() ,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")


SampleNum <- length(dataRAW) - 1
data <- dataRAW #[,c(2,13:length(dataRAW))]
names(data) <- c("Chromosome","ESC 2i","EpiLC","PGCLC","Epi E6.5","PGC E9.5","PGC E10.5", "PGC E11.5", "PGC E13.5f", "PGC E13.5m", "ESC 2i", "ESC ser")
#data <- data[,c(1,2,3,5,4,6,7,8)]

#data <- data[1:100,]

#Exclude non-quantified samples
data <- data[apply(data[,2:(1+SampleNum)],1,min)>=0,]
data_matrix <-data.matrix(data[,2:(1+SampleNum)])

data_mean <- apply(data_matrix,2,mean)
data_median <- apply(data_matrix,2,median)
data_sd <- apply(data_matrix,2,sd)

c("white","grey","PGCLC","Epi E6.5","PGC E9.5","PGC E10.5", "PGC E11.5", "PGC E13.5f", "PGC E13.5m", "ESC 2i", "ESC ser")

colUsed <- c("white","grey","#fe1a00","#243e3e","#fed10a","#fed10a","#fed10a","#fed10a","#fed10a","white","grey")

boxplot(data_matrix,
        axes=FALSE,axisnames=FALSE,
        col=colUsed,
        main="CpG methylation",ylab="% CpG methylation")
axis(1, at=c(1:SampleNum),labels = colnames(data_matrix),tick=T)
axis(2, at=seq(0 , 100, 10))

library(vioplot)
vioplot(data[,2],data[,3],data[,4],data[,5],data[,6],data[,7],data[,8],
        col="gold")
axis(1, at=c(1:SampleNum),labels = colnames(data_matrix),tick=T)
axis(2, at=seq(0 , 100, 10),ylab="% CpG methylation")




bp <- barplot(data_mean,
              axes=FALSE,axisnames=TRUE,ylim=c(0,100),
              col=colUsed,
              main="Average CpG methylation",ylab="% CpG methylation")

# The x-axis
#axis(1, at = bp,labels = F,line=NA)
# The y-axis 
axis(2, at=seq(0 , 100, 10))

# Now plot the error bars
# 1. The lower bar &  the vertical bar
#segments(bp, data_mean, bp, data_mean - data_sd, lwd=2)
#segments(bp - 0.1, data_mean - data_sd, bp + 0.1, data_mean - data_sd, lwd=2)

# 2. The upper bar & the vertical bar
segments(bp, data_mean, bp, data_mean + data_sd, lwd=2)
segments(bp - 0.1, data_mean + data_sd, bp + 0.1, data_mean + data_sd, lwd=2)



