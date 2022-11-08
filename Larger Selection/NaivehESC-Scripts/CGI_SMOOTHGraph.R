require(grid)

setwd("~/Documents/Cambridge/Projects/huESC Ge et al/2015Oct/CGI/")
filenameRAW <- "CGI5_5Methylation.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
dataRAW <- dataRAW[,c(2,13,14,15,16,17)]
names(dataRAW) <- c("Chromosome","ICM (Guo)","HNES1","Reset H9","HNES1 primed","H9")

SampleNum <- 5

chrom <- "All"

N <- dataRAW
N <- dataRAW[dataRAW$Chromosome==chrom,]

#TESTING
i<-3
j<-5
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
      
      #Create PDF with Frame
      pdf(paste("CGI5_5_Chr",chrom,"_",names(N)[i],"_vs_",names(N)[j],".pdf",sep=""),width=5,height=5,useDingbats=T,colormodel="cmyk")
      par(pty="s")
      smoothScatter(x2~x1,
                     bandwidth=4, nrpoints=0,
                     xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,
                     xlab = names(N)[i],ylab=names(N)[j],main=paste("%CGI methylation\n","Chromosome",chrom)
      )
      plot(x2~x1,
                     pch=19, cex=.5,type ="p",col="blue",
                     xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,
                     xlab = names(N)[i],ylab=names(N)[j],main=paste("%CGI methylation\n","Chromosome",chrom)
      )
 
      dev.off()
    }
  }
}