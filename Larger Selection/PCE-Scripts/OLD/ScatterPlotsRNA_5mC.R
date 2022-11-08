#library(reshape2)
#library(plyr)
#require(grid)
colorplot <- function (data.plot, title="Methylation vs Expression",BSLab = "%CpG methylation",PDF=FALSE,Reg=TRUE)  {
  # Remove incomplete rows & Remove non expressed genes & 0% methylation 
  
  data.plot <- data.plot[data.plot[2]>min(data.plot[2]),]
  data.plot <- data.plot[data.plot[3]!=0,]
  names(data.plot) <- c("Gene","RNA","BS")
  
  #Making plots ...
  RSQ<- round(cor.test(data.plot$BS, data.plot$RNA)[4][[1]],3)
  x1 <- data.plot$BS
  x2 <- data.plot$RNA
  df <- data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  
  #Create PDF 
  if (PDF) {pdf(paste("RNAvsBS_",title,".pdf",sep=""),width=5,height=5,useDingbats=FALSE)}
  #Plot
  plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
       axes = T,frame.plot = F,
       xlim = c(0,100),
       xlab = BSLab ,ylab="rel. mRNA expression",
       main=title)
  if (Reg) {
    abline(lm(x2~x1), col="black") # regression line (y~x) 
    text(85,max(x2)-2,paste("R =",round(RSQ,3)))
  }
  if (PDF) {dev.off()}
}

setwd("~/Documents/Cambridge/DataSeq/201510_Islets/")
filenameRAWRNA <- "RNASeq_Islets/RNASeq_Islets_OBExpValues.txt"
  dataRAWRNA <- read.csv(filenameRAWRNA,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
  dataRAWRNA <- cbind(Gene=dataRAWRNA$Probe,dataRAWRNA[,c(13,14,15)])
  names(dataRAWRNA) <- c("Gene","OB_RNA","OK_RNA","WT_RNA")
  head(dataRAWRNA)
  
  
filenameRAWBS <- "BSSeq_Islets/Promoter_CGI_Methylation_OB.txt"
  dataRAWBS <- read.csv(filenameRAWBS,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
  dataRAWBS <- cbind(Gene=sub("_upstream","",as.character(dataRAWBS$Probe)),dataRAWBS[,c(13,14,15)],CGI=(dataRAWBS$Feature != "null"))
  #dataRAWBS <- dataRAWBS[,c(1,13,14,15)]
  names(dataRAWBS) <- c("Gene","OB_BS","OK_BS","WT_BS","CGI")
  head(dataRAWBS)
  

#dataRAWRNA <- dataRAWRNA[order(dataRAWRNA[,1]),]
#rownames(dataRAWRNA) <- c(1:length(dataRAWRNA[,1]))
#dataRAWBS <- dataRAWBS[order(dataRAWBS[,1]),]
#rownames(dataRAWBS) <- c(1:length(dataRAWBS[,1]))

#dataRAWRNA<-dataRAWRNA[c(1:100),]
#dataRAWBS<-dataRAWBS[c(1:100),]

data <- merge(dataRAWRNA,dataRAWBS,by.x='Gene',by.y="Gene")
head(data)
#data <- na.omit(data)
dataCGI <- data[data$CGI,]
dataNONCGI <- data[!data$CGI,]
head(dataNONCGI)
head(dataCGI)

#Split into comparsions

#OB/OB
  colorplot(na.omit(dataNONCGI[,c(1,2,5)]),"OBOB Genebody",BSLab = "Genebody methylation", PDF=F,Reg=T)
#OB Keto
   colorplot(na.omit(dataCGI[,c(1,3,6)]),"OB Keto Genebody",BSLab = "Genebody methylation", PDF=F,Reg=T)
#WT
  colorplot(na.omit(dataCGI[,c(1,4,7)]),"WT Genebody",BSLab = "Genebody methylation", PDF=F,Reg=T)
  
  
  
  
  


