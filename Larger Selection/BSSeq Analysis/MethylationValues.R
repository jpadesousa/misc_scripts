setwd("~/Documents/Cambridge/DataSeq/mESC_Ser2i/")
filenameRAW <- file.choose() 
dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
SampleNum <- length(dataRAW) - 12
dataRESC4 <- dataRAW[,c(13:length(dataRAW))]
#names(data) <- c("Serum","Serum","2i 24h","2i 24h","2i LT","2i LT")

data <- dataBG
names(data) <- c("Serum - Random Probe","2i - Random Probe","Serum - Low H3K9me2","2i - Low H3K9me2","Serum - High H3K9me2","2i - High H3K9me2")

#Exclude non-quantified samples
#data <- data[apply(data[,2:(1+SampleNum)],1,min)>=0,]
data_matrix <-data.matrix(data)


data_mean <- apply(data_matrix,2,mean)
data_median <- apply(data_matrix,2,median)
data_sd <- apply(data_matrix,2,sd)


colUsed <- c("darkblue","darkblue","blue","blue")

boxplot(dataBG[,1],dataRESC4[,1],dataBG[,2],dataRESC4[,2],
        axes=F,axisnames=F,
        col=colUsed,notch=F,
        main="H3K9me3",ylab="H3K9me3 enrichment")
axis(1, at=c(1:4),labels = c("Random","Res","Random","Res"),tick=T)
axis(2, at=seq(-10 , 15, 5))

boxplot(dataBG[,3],dataRESC4[,3],dataBG[,4],dataRESC4[,4],
        axes=F,axisnames=F,
        col=colUsed,notch=F,
        main="H3K9me2",ylab="H3K9me2 enrichment")
axis(1, at=c(1:4),labels = c("Random","Res","Random","Res"),tick=T)
axis(2, at=seq(-10 , 15, 5))



t.test(dataBG[,3],dataRESC4[,3])
serTT$p.value
TT2i24h <- t.test(data[,3],data[,4])
TT2i24h$p.value
TT2iLT <- t.test(data[,5],data[,6])
TT2iLT$p.value

library(beanplot)
beanplot(data[,],
        col="gold")
axis(1, at=c(1:SampleNum),labels = colnames(data_matrix),tick=T)
axis(2, at=seq(0 , 100, 10),ylab="% CpG methylation")




bp <- barplot(data_mean,
              axes=FALSE,axisnames=TRUE,ylim=c(0,100),
              col=colUsed,
              main="Average CpG methylation",ylab="% CpG methylation")

# The x-axis

axis(1, at = bp,labels = data2_names,line=NA,tick = F)
# The y-axis 
axis(2, at=seq(0 , 100, 10))

# Now plot the error bars
# 1. The lower bar &  the vertical bar
#segments(bp, data_mean, bp, data_mean - data_sd, lwd=2)
#segments(bp - 0.1, data_mean - data_sd, bp + 0.1, data_mean - data_sd, lwd=2)

# 2. The upper bar & the vertical bar
segments(bp, data_mean, bp, data_mean + data_sd, lwd=2)
segments(bp - 0.1, data_mean + data_sd, bp + 0.1, data_mean + data_sd, lwd=2)



