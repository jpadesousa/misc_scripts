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
color.plot <-   function(MethVal,HistVal,HistName,plotting="r",plotName="test",cex=0.5){
  x1 <- MethVal[[1]]
  x2 <- MethVal[[2]]
  df <- data.frame(x1,x2)
  
  # Colors by Histones    
  df$dens <- round(rescale(HistVal,to = c(1, 256)))
  
  ## Map to colors
  cols <-  colorRampPalette(c("#000099","#00FEFF","#45FE4F","#FCFF00","#FF9400","#FF3100"))(256)
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
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,bty ="n",
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
color.hist <-   function(MethVal,HistVal,HistName,plotting="r",plotName="test",cex=0.5){
  x1 <- MethVal[[1]]
  x2 <- MethVal[[2]]
  df <- data.frame(x1,x2)
  
  # Colors by Histones    
  df$dens <- round(rescale(HistVal,to = c(1, 256)))
  
  ## Map to colors
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
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,bty ="n",
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
color.plot2C <- function(MethVal,HistVal,HistName,plotting="r",plotName="test",cex=0.5){
  x1 <- MethVal[[1]]
  x2 <- MethVal[[2]]
  df <- data.frame(x1,x2)
  
  # Colors by Histones    
  df$dens <- round(rescale(HistVal,to = c(1, 3)))
  
  ## Map to colors
  cols <-  colorRampPalette(c("#DDDDDD","dark blue","#FF3100"))(3)
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
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,bty ="n",
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
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,bty ="n",
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
density.hist <- function(MethVal,HistBool,HistName,plotting="r",plotName="test",cex=0.5){
  n1 <- MethVal[[1]][HistBool == 0]
  n2 <- MethVal[[2]][HistBool == 0]
  
  # RSQ <- round(cor.test(x1, x2)[4][[1]],3)  #Calculate RSQ of comparison
  
  # Colors by density of Mark present
  ## Use densCols() output to get density at each point
  d1 <- MethVal[[1]][HistBool == 1]
  d2 <- MethVal[[2]][HistBool == 1]
  d <- densCols(d1,d2, colramp=colorRampPalette(c("black", "white")))
  
  df <- data.frame(d1,d2)
  df$dens <- col2rgb(d)[1,] + 1L
  
  ## Map to colors
  cols <-  colorRampPalette(c("#000099","#00FEFF","#45FE4F","#FCFF00","#FF9400","#FF3100"))(256)
  df$col <- cols[df$dens]
  
  ## Create plot
  if (plotting == "r") {
    par(pty="s")
    plot(n2~n1, pch=19, col="#DDDDDD", cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,bty ="n",
         xlab = names(MethVal)[1],ylab=names(MethVal)[2],
         main=paste("Density of ",HistName))
    par(new=TRUE)
    plot(d2~d1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,ann=F)
  }
  
  ## Create PDF
  if (plotting == "p") {
    pdf (paste(plotName,"_",HistName,".pdf", sep=""), width = 10, height = 10, useDingbats=FALSE,colormodel='rgb')
    par(pty="s")
    plot(n2~n1, pch=19, col="#DDDDDD", cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,bty ="n",
         xlab = names(MethVal)[1],ylab=names(MethVal)[2],
         main=paste("Density of ",HistName))
    par(new=TRUE)
    plot(d2~d1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,ann=F)
    dev.off()
  }
  
  ## Create TIFF
  if (plotting == "t") {
    tiff(paste(plotName,"_",HistName,".tiff", sep=""), width = 10, height = 10, units = 'in', res = 300, compression = 'lzw')
    par(mar=c(0,0,0,0),mai=c(0,0,0,0), xpd = NA,pty="s") 
    plot(n2~n1, pch=19, col="#DDDDDD", cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,bty ="n",
         xaxs = "i", yaxs = "i")
    par(new=TRUE)
    plot(d2~d1, data=df[order(df$dens),], pch=19, col=col, cex=cex,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,bty ="n",
         xaxs = "i", yaxs = "i")
    dev.off()
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
onoff.histones <- function(HistVal, ByPer = T, Value=50){
  #Divide Histone Values
  NewHistVal <- HistVal
  for (i  in 1:ncol(HistVal)) {
    if (ByPer) {Seperator <- quantile(t(HistVal[,i]),Value/100)}
    else {Seperator <- Value}
    NewHistVal[(HistVal[,i] <= Seperator),i] <- 0
    NewHistVal[(HistVal[,i] > Seperator),i] <- 1
  }
  return (NewHistVal)
}
########

########


setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/humBSSeq/Scatter/")
filenameRAW <- "humBSData_50CpG_5_1.txt"
testRead <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",nrows = 5)
#classes <- sapply(testRead, class)
classes <- c("character","factor",rep("NULL",10),rep("numeric",ncol(testRead)-12))
MethValues <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",colClasses = classes)

Samples <- ncol(MethValues) - 2

filenames <- list.files(path = "Annotation", pattern="*.txt", full.names=F)
for (filenameRAW in filenames) {
  #filenameRAW <- "HighD12.txt"
  classes <- c("character","factor","integer","integer","numeric")
  Input <- read.csv(paste("Annotation/",filenameRAW,sep=""),sep="\t",
                    header=TRUE,skip=0,fill=TRUE,dec = ".",colClasses = classes)

  Input[,1] <- paste("Chr",Input[,2],":",Input[,3],"-",Input[,4],sep="")
  match(Input[,1],MethValues[,1]) -> pos

  MethValues <- cbind(MethValues,rep(0,nrow(MethValues)))
  names(MethValues)[ncol(MethValues)] <- strsplit(filenameRAW,".txt")
  MethValues[pos,ncol(MethValues)] <- 1
  }

#head(MethValues,n=10)

#Combine Highd12 & LowD12
MethValues$HighLow <- MethValues$HighD12 + 0.5 * MethValues$LowD12

Marks <- ncol(MethValues) - Samples - 2

#head(MethValues[,c(9,13)])
#head(MethValues[,c((Samples+3):(Samples+Marks+2))])

## START FROM HERE
N <- na.omit(MethValues[,c(9,13)])                                            #Methylation comparison
M <- MethValues[as.integer(rownames(N)),c((Samples+3):(Samples+Marks+2))]     #Overlay marks
#Scale Histone Values
# M <- scale.histones(M,50,99)


pos <- sample(nrow(N),10000) #TESTING
cex = 0.5
pos <- rownames(N) # Use all positions
cex = 0.1
mark = 4

names(M)[mark]
color.plot(N[pos,],M[pos,mark],names(M)[mark],cex=cex)    #Label by Mark
color.hist(N[pos,],M[pos,mark],names(M)[mark],cex=cex)    #Label by Mark, low values in grey
color.plot2C(N[pos,],M[pos,mark],names(M)[mark],cex=cex)  #Label by Mark with 3 Values, ie NO - Mark1 - Mark2
density.hist(N[pos,],M[pos,mark],names(M)[mark],cex=cex)  #Label by Density of Mark

density.plot(N[pos,],cex=cex)                             #Label by Density


plotName="humBSData_50CpG_5_1"


density.hist(N[pos,],M[pos,mark],names(M)[mark],cex=cex,plotting="t",plotName=plotName)  #Label by Density of Mark


density.plot(N[pos,],cex=0.5,plotting="p",plotName=plotName)
density.plot(N,cex=0.1,plotting="t",plotName=plotName)
#color.bar (colorRampPalette(c("#000099","#00FEFF","#45FE4F","#FCFF00","#FF9400","#FF3100"))(256),0,1,6)


for (mark in 1:2)  color.plot2C(N[pos,],M[pos,mark],names(M)[mark],cex=0.5,plotting="p",plotName=plotName)
for (mark in 1:2) color.plot2C(N[pos,],M[pos,mark],names(M)[mark],cex=0.1,plotting="t",plotName=plotName)

#color.bar (colorRampPalette(c("#DDDDDD","#DDDDDD","#DDDDDD","#FCFF00","#FF9400","#FF3100"))(256),0,1,6)
                                 


M <- Input[as.integer(rownames(N)),]     #Overlay marks
low50 <- quantile(t(M),1/100)
top1 <- quantile(t(M),99/100)
M[(M < low50)] <- low50
M[(M > top1)] <- top1

pos <- sample(nrow(N),100000) #TESTING
pos <- rownames(N) # Use all positions

color.plot(N[pos,],M[pos],"CG density",cex=cex,plotName=plotName,plotting="t")    #Label by Mark
