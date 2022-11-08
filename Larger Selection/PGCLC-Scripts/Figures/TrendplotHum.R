setwd("/Users/meyennf/Desktop/Trendplot/")


filenameRAW <- "humGene_2.txt"
Input <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")


Pos <- c(1:5000,seq(from = 5050, to = 20000, length.out = 200),20001:25000)
cols <-  colorRampPalette(c("#000000","#00FEFF","#45FE4F","#FCFF00","#FF9400","#FF3100"))(ncol(Input))


pdf ("Human Gene TrendPlot2.pdf", width = 10, height = 10, useDingbats=FALSE,colormodel='rgb')

par(cex.lab=1, cex.axis=0.7, cex.sub=0.5, mar=c(5,5,3,3))

plot(Pos,rep(0,length(Pos)),type="n",
     main = "Human Gene TrendPlot",
     ylab = "CpG methylation (%)", xlab ="",
     ylim = c(0,80)
     )

for (i in 2:ncol(Input)) {
  par(new=TRUE)
  plot(Pos,Input[,i],
       type="l", 
       ylim = c(0,80),ann=FALSE, axes=FALSE,col=cols[i-1])
}
legend("bottomright", legend = names(Input)[-1], 
       col = cols, lty= 1, lwd = 5, cex=0.5)          


dev.off()
