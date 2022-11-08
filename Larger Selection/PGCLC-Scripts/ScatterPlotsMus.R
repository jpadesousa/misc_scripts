require(grid)
require(scales)

######
color.bar <- function(lut, min=-1, max=-min, nticks=11,ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
########
color.plot <- function(MethVal,HistVal,HistName,plotting="r",plotName="test",cex=0.5){
  x1 <- MethVal[[1]]
  x2 <- MethVal[[2]]
  df <- data.frame(x1,x2)
  
  # Colors by Histones    
  df$dens <- round(rescale(HistVal,to = c(1, 256)))
  
  ## Map to colors
  #cols <-  colorRampPalette(c("#000099","#00FEFF","#45FE4F","#FCFF00","#FF9400","#FF3100"))(256)
  cols <-  colorRampPalette(c("#DDDDDD","#DDDDDD","#DDDDDD","#FCFF00","#FF9400","#FF3100"))(256)
  df$col <- cols[df$dens]
  
  ## Create PDF
  if (plotting == "p") {
    pdf (paste(plotName,"_",HistName,".pdf", sep=""), width = 10, height = 10, useDingbats=FALSE,colormodel='rgb')
    par(pty="s")
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,bty ="n",
         xlab = names(MethVal)[1],ylab=names(MethVal)[2],
         main=paste("Colored by",HistName))
    dev.off()
  }
  
  ## Create TIFF
  if (plotting == "t") {
    tiff(paste(plotName,"_",HistName,".tiff", sep=""), width = 10, height = 10, units = 'in', res = 300, compression = 'lzw')
    par(mar=c(0,0,0,0),mai=c(0,0,0,0), xpd = NA,pty="s") 
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.5,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,#bty ="n",
         xaxs = "i", yaxs = "i")
    dev.off()
  }
  
  ## Create plot
  if (plotting == "r"){
    par(pty="s")
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
       xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,bty ="n",
       xlab = names(MethVal)[1],ylab=names(MethVal)[2],
       main=paste("Colored by",HistName))
  }
}
########
density.plot <- function(MethVal,plotting="r",plotName="test",cex=0.5){
  x1 <- MethVal[[1]]
  x2 <- MethVal[[2]]
  df <- data.frame(x1,x2)
  
  # RSQ <- round(cor.test(x1, x2)[4][[1]],3)  #Calculate RSQ of comparison
  
  
  # Colors by density    
  ## Use densCols() output to get density at each point
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map to colors
  cols <-  colorRampPalette(c("#000099","#00FEFF","#45FE4F","#FCFF00","#FF9400","#FF3100"))(256)
  df$col <- cols[df$dens]
  
  ## Create PDF
  if (plotting == "p") {
    pdf (paste(plotName,"_Scatter.pdf", sep=""), width = 10, height = 10, useDingbats=FALSE,colormodel='rgb')
    par(pty="s")
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,bty ="n",
         xlab = names(MethVal)[1],ylab=names(MethVal)[2],
         main=paste("%CpG methylation"))
    dev.off()
  }

  ## Create TIFF
  if (plotting == "t") {
    tiff(paste(plotName,"_Scatter.tiff",sep=""), width = 10, height = 10, units = 'in', res = 300, compression = 'lzw')
    par(mar=c(0,0,0,0),mai=c(0,0,0,0), xpd = NA,pty="s") 
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.5,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,#bty ="n",
         xaxs = "i", yaxs = "i")
    
    dev.off()
  }
  
  ## Create plot
  if (plotting == "r") {
    par(pty="s")
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
       xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,bty ="n",
       xlab = names(MethVal)[1],ylab=names(MethVal)[2],
       main=paste("%CpG methylation"))
  }

}
########
scale.histones <- function(HistVal,lowest=50,highest=99){
  #Scale Histone Values
  for (i  in 1:ncol(HistVal)) {
    #    hist(HistVal[,i], main=paste(names(HistVal)[i],"\nRange: ",min(HistVal[,i])," <> ",max(HistVal[,i]),"\nMean: ",apply(t(HistVal[,i]),1,mean)))
    low50 <- quantile(t(HistVal[,i]),lowest/100)
    top1 <- quantile(t(HistVal[,i]),highest/100)
    HistVal[(HistVal[,i] < low50),i] <- low50
    HistVal[(HistVal[,i] > top1),i] <- top1
    #    hist(HistVal[,i], main=paste(names(HistVal)[i],"\nRange: ",min(HistVal[,i])," <> ",max(HistVal[,i]),"\nMean: ",apply(t(HistVal[,i]),1,mean)))
  }
  return (HistVal)
}
########


setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musBSSeq/ChIPSeq/")

filenameRAW <- "HistoneValues_Mus.txt"
HistValues <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
HistValues <- HistValues[,c(1,2,13:ncol(HistValues))]

filenameRAW <- "HistoneValues_Mus_H3K9me3.txt"
Input <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
Input <- Input[,c(13:ncol(Input))]
HistValues <- cbind(HistValues,Input)
remove(Input)

names(HistValues) <- c("Probe","Chromosome","H3K27ac","H3K4me3","H3K27me3","H3K9me2","pan H3","H3K9me3")

filenameRAW <- "musBSData_50CpG_5_1.txt"
MethValues <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
MethValues <- MethValues[,c(1,2,13:ncol(MethValues))]


IAPs <- read.csv("IAPOverlap.txt",sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
IAPs <- IAPs[,c(1,2,7)]
cbind(IAPs,IAP=rep(0,nrow(IAPs))) -> IAPs
IAPs[IAPs[,3] != "null",4] <- 1

HistValues <- cbind(HistValues,IAP=IAPs[,4])
remove(IAPs)

## START FROM HERE
N <- na.omit(MethValues[,c(8,10)])                                #Methylation comparison
H <- HistValues[as.integer(rownames(N)),c(3:ncol(HistValues))]    #Histone marks

#Scale Histone Values
H <- scale.histones(H,50,99)

#TESTING
pos <- sample(nrow(N),50000)

color.plot(N[pos,],H[pos,7],names(H)[7],cex=1)

density.plot(N,cex=0.1,plotting="r",plotName="musBSData_50CpG_5_1")
density.plot(N[pos,],cex=0.5,plotting="p",plotName="musBSData_50CpG_5_1")
density.plot(N,cex=0.1,plotting="t",plotName="musBSData_50CpG_5_1_all")

#for (Histmark in 1:6)  color.plot(N[pos,],H[pos,Histmark],names(H)[Histmark],cex=1,plotting="p",plotName="musBSData_50CpG_5_1")
for (Histmark in 1:7)  color.plot(N[pos,],H[pos,Histmark],names(H)[Histmark],cex=1,plotting="t",plotName="musBSData_50CpG_5_1")

for (Histmark in 1:7)  color.plot(N,H[,Histmark],names(H)[Histmark],cex=0.1,plotting="t",plotName="musBSData_50CpG_5_1_01")

color.bar (colorRampPalette(c("#DDDDDD","#DDDDDD","#DDDDDD","#FCFF00","#FF9400","#FF3100"))(256),0,1,6)


