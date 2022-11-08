setwd("/Users/meyennf/Desktop/PGC/humData//")
library(beanplot)

column.count = max(count.fields("AllProbes.txt", sep = "\t"))
AllProbes <- read.delim("AllProbes.txt",fill=F,header=T,stringsAsFactors = FALSE,colClasses = c(rep("NULL",12),rep(NA,column.count-12)))
LowProbes <- read.delim("LowProbes.txt",fill=F,header=T,stringsAsFactors = FALSE,colClasses = c(rep("NULL",12),rep(NA,column.count-12)))
HighProbes <- read.delim("HighProbes.txt",fill=F,header=T,stringsAsFactors = FALSE,colClasses = c(rep("NULL",12),rep(NA,column.count-12)))
i=1
for (i in 1:(column.count-12) )  
  {
  pdf(paste(names(AllProbes[i]),".pdf",sep=""),width=4)
  par(las=1,mar=c(1,5,3,1),cex.axis=0.8)
  beanplot(c(AllProbes[i],LowProbes[i],HighProbes[i]),
           bw="nrd0",
           ylim=c(0,100),
           what=c(0,1,1,0),
           show.names=F,
           ylab="Methylation percentage")
  title(names(AllProbes[i]))
  dev.off()
  }

plot.amplicon <- function(condition.name) {
  
  data.subset <- data[data[[2]]==condition.name,]
  
  # Remove lines with no data
  data.subset[apply(data.subset[,3:ncol(data.subset)],1,function(y) return(sum(!is.na(y))))>10,] -> data.subset
  
  if (nrow(data.subset) == 0) {
    return(0)
  }

  
  
  
  row.names(data.subset) <- data.subset[[1]]
  
  par(las=1,mar=c(5,8,4,2),cex.axis=0.8)
  beanplot(as.data.frame(t(data.subset[,3:ncol(data.subset)])),
           bw="nrd0",
           horizontal = TRUE,
           ylim=c(0,100),
           what=c(0,1,0,0),
           xlab="Methylation percentage")
  title(paste("Methylation distributions for",condition.name))
  
}

sapply(unique(data[[2]]),plot.amplicon)


