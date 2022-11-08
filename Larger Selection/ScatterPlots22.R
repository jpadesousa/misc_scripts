require(grid)

setwd("/Users/meyennf/Documents/Cambridge/Projects/hESC/Guo_2/")  ## Data folder
filenameRAW <- "500kb_5mC.txt"                                    ## Data file

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
dataRAW <- dataRAW[,c(2,13:length(dataRAW))]

## Change names if required

#names(dataRAW) <- c("Chromosome","ICM Guo et al. 2014","HNES1 Guo et al. 2016","Shef6-EOS (primed)","Shef6-EOS (CR p5)","Shef6-EOS (CR p9)","Shef6-EOS (CR LT)","H9-NK2 (primed) Takashima et al. 2014","H9-NK2 (reset) Takashima et al. 2014")

N <- dataRAW

## For TESTING purpose: select only one comparison
i<-4
j<-3

## For TESTING purpose: reduce the amount of data
#N <- dataRAW[1:100,]

## Inititalize drawing system
dev.off()
plot(1)

#Generate plots for all possible combinations
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
      cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
      df$col <- cols[df$dens]
    
      ## Create plot, reordering rows so that densest points are plotted on top

      ## Create PDF with Frame in DataFolder - Change to png or tiff or pdf depending or requirements
      pdf(paste("5mC_Scatter_",names(N)[i],"_vs_",names(N)[j],".pdf",sep=""),width=5,height=5,useDingbats=FALSE,colormodel="rgb")
      par(pty="s")
      plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
           xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,#bty ="n",
           xlab = names(N)[i],ylab=names(N)[j],
           main="%CpG methylation")
   
      ## Add Correlation and line of ideal correlation   
      #abline(c(1,1), col="black") 
      #abline(lm(x2~x1), col="black") # regression line (y~x)
      #text(95,0,paste("R =",round(RSQ,3)),cex=0.4)
      dev.off()
      
    }
  }
}