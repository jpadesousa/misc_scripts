require(grid)

setwd("~/Documents/Cambridge/DataSeq/201510_Islets/RNASeq_Islets")
filenameRAW <- "RNASeq_Islets_OBExpValues.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
dataRAW <- dataRAW[,c(1,13,14,15)]
names(dataRAW) <- c("Gene","OBOB","OBOB Keto","wildtype")

SampleNum <- 3


N <- dataRAW

#TESTING
#i<-2
#j<-3
#data <- data[1:100,]
#TESTING

#Inititalize drawing system
#dev.off()
#plot(1)

for (i in 2:length(N)){
  for (j in 2:length(N)){
    x1 <- N[[i]]
    x2 <- N[[j]]
    RSQ <- round(cor.test(x1, x2)[4][[1]],3)
    if (RSQ != 1) {
      df <- data.frame(x1,x2)
    
      ## Use densCols() output to get density at each point
      x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
      df$dens <- col2rgb(x)[1,] + 1L
      
      ## Map densities to colors
      cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                                  "#FCFF00", "#FF9400", "#FF3100"))(256)
      df$col <- cols[df$dens]
      
      xmin <- min(x1,x2)
      xmax <- max(x1,x2)
      
      ## Create plot, reordering rows so that densest points are plotted on top
#       dev.new()
#       #op <- par(mar = rep(0, 4))
#       par(mar=c(0,0,0,0),mai=c(0,0,0,0),pty="s",cex=0.8) 
#       plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
#            xlim=c(xmin,xmax),ylim=c(xmin,xmax),
#            axes = F,frame.plot = F,asp=1,#bty ="n",
#            xaxs = "i", yaxs = "i" #,xlab = "",ylab=""
#            )
#       capture <- grid.cap()
#       dev.off()
    
      #Create PDF with plot & Frame
      pdf(paste("RNASeq_",names(N)[i],"_vs_",names(N)[j],".pdf",sep=""),width=5,height=5,useDingbats=FALSE)
      par(pty="s",cex=0.8)
      plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
           xlim=c(xmin,xmax),ylim=c(xmin,xmax),
           axes = T,frame.plot = F,asp=1,
           xlab = names(N)[i],ylab=names(N)[j])
      # abline(c(1,1), col="black") 
      abline(lm(x2~x1), col="black") # regression line (y~x)
      text(xmax-4,xmin,paste("R =",round(RSQ,3)))
      dev.off()
      
    }
  }
}
