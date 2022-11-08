require(grid)

setwd("~/Documents/Cambridge/Projects/huESC Ge et al/2015Oct/Data/")
filenameRAW <- "500kbCpG_1_1_BSPipeline.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
dataRAW <- dataRAW[,c(2,13,14,15,16,17,18,19)]
names(dataRAW) <- c("Chromosome","ICM (Guo)","HNES1","HNES3","Reset H9","HNES1 primed","HNES3 primed","H9")

SampleNum <- 7


chrom <- "6"
chromList <- levels(dataRAW$Chromosome)

dataNAIVE <- dataRAW[,c(1,2,3)]

N <- dataRAW[dataNAIVE$Chromosome==chrom,]

#TESTING
i<-3
j<-2
#N <- dataRAW[1:100,]
#TESTING

#Inititalize drawing system
dev.off()
plot(1)

for (chrom in chromList){
    N <- dataRAW[dataNAIVE$Chromosome==chrom,]
    x1 <- N[[i]]
    x2 <- N[[j]]
    RSQ <- round(cor.test(x1, x2)[4][[1]],3)
    df <- data.frame(x1,x2)
      
      ## Use densCols() output to get density at each point
      x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
      df$dens <- col2rgb(x)[1,] + 1L
      
      ## Map densities to colors
      cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                                  "#FCFF00", "#FF9400", "#FF3100"))(256)
      df$col <- cols[df$dens]
      
      
      ## Create plot, reordering rows so that densest points are plotted on top
      dev.new(units="px", width=1600, height=1600, res=300)
      #op <- par(mar = rep(0, 4))
      par(mar=c(0,0,0,0),mai=c(0,0,0,0), xpd = NA,pty="s") 
      plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
           xlim=c(0,100),ylim=c(0,100),axes = F,frame.plot = F,asp=1,#bty ="n",
           xaxs = "i", yaxs = "i" #,xlab = "",ylab=""
      )
      capture <- grid.cap()
      dev.off()
      
      #Create PDF with Frame
      pdf(paste("CpG_Chr",chrom,"_",names(N)[i],"_vs_",names(N)[j],"RSQ.pdf",sep=""),width=5,height=5,useDingbats=FALSE,colormodel="cmyk")
      par(pty="s")
      plot(1,type="n",
           xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,
           xlab = names(N)[i],ylab=names(N)[j],main=paste("%CpG methylation\n","Chromosome",chrom))
      rasterImage(capture,0,0,100,100)
      #abline(c(1,1), col="black") 
      #abline(lm(x2~x1), col="black") # regression line (y~x)
      text(95,0,paste("R =",round(RSQ,3)),cex=0.4)
      dev.off()
      
    }
  }
}