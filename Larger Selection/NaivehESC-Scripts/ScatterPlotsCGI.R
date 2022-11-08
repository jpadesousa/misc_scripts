require(grid)

setwd("~/Documents/Cambridge/Projects/huESC/Guo_2/")
filenameRAW <- "CGI_5mC.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
dataRAW <- dataRAW[,c(2,13:21)]
names(dataRAW) <- c("Chromosome","ICM (Guo)",
                    "S6_EOS_KSR","S6_EOS_p5","S6_EOS_p9","S6_EOS_p15",
                    "HNES1","HNES3","Conventional H9","naive H9")

N <- dataRAW

#TESTING
i<-4
j<-3
#N <- dataRAW[1:100,]
#TESTING

#Inititalize drawing system
dev.off()
plot(1)

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

      
      ## Create plot, reordering rows so that densest points are plotted on top

      #Create PDF with Frame
      pdf(paste("CGI//CGI5mC_Scatter_",names(N)[i],"_vs_",names(N)[j],".pdf",sep=""),width=5,height=5,useDingbats=FALSE,colormodel="rgb")
      par(pty="s")
      plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
           xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,#bty ="n",
           xlab = names(N)[i],ylab=names(N)[j],main=paste("CGI methylation")
           )
      
      abline(c(1,1), col="black") 
      #abline(lm(x2~x1), col="black") # regression line (y~x)
      text(95,0,paste("R =",round(RSQ,3)),cex=0.4)
      dev.off()
      
    }
  }
}