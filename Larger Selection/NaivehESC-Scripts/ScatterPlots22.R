require(grid)

setwd("/Users/meyennf/Documents/Cambridge/Projects/hESC/Guo_2/_FinalFigures/")
filenameRAW <- "Sel_Val_CGIProm.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
dataRAW <- dataRAW[,c(2,15:length(dataRAW))]
names(dataRAW) <- c("Chromosome","ICM Guo et al. 2014","HNES1 Guo et al. 2016","HNES1 (primed) Guo et al. 2016","Shef6-EOS (primed)","Shef6-EOS (CR p5)","Shef6-EOS (CR p9)","Shef6-EOS (CR LT)","H9-NK2 (primed) Takashima et al. 2014","H9-NK2 (reset) Takashima et al. 2014")
#names(dataRAW) <- c("Chromosome","ICM Guo et al. 2014","HNES1 Guo et al. 2016","Shef6-EOS (primed)","Shef6-EOS (CR p5)","Shef6-EOS (CR p9)","Shef6-EOS (CR LT)","H9-NK2 (primed) Takashima et al. 2014","H9-NK2 (reset) Takashima et al. 2014")

names(dataRAW)
chrom <- "Sel"

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
    
    #  df$col <- "red" 
      
      
      ## Create plot, reordering rows so that densest points are plotted on top

      #Create PDF with Frame
      pdf(paste("Scatter-CGIPromSel/Sel_5mC_Scatter_",names(N)[i],"_vs_",names(N)[j],".pdf",sep=""),width=5,height=5,useDingbats=FALSE,colormodel="rgb")
      par(pty="s")
      plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
           xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,#bty ="n",
           xlab = names(N)[i],ylab=names(N)[j],
           main=paste("%CpG methylation\n","Chromosome",chrom)
           #main=paste("CGI Promoter methylation")
           )
      
      #abline(c(1,1), col="black") 
      #abline(lm(x2~x1), col="black") # regression line (y~x)
      #text(95,0,paste("R =",round(RSQ,3)),cex=0.4)
      dev.off()
      
    }
  }
}