setwd("/Users/meyennf/Desktop/Imprint/Imprint")
library(beanplot)

column.count = max(count.fields("amplicon_methylation_percentages.txt", sep = "\t"))

read.delim("amplicon_methylation_percentages.txt",col.names=1:column.count,fill=TRUE,header=FALSE,stringsAsFactors = FALSE) -> data
gsub("lane\\d+_[GATC]+_","",data[[1]]) -> data[[1]]
gsub("_L\\d+_R\\d+_val_\\d+.fq.gz_bismark_bt2_pe_merged.bam","",data[[1]]) -> data[[1]]

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
