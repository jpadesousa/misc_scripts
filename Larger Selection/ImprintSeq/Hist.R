impHist <- function (x) {
  x <- unlist(x)
  sample <- substr(x[1],19,nchar(x[1])-46)
  gene <- x[2]
  if( length(x)>2) {
    val <- as.numeric(x[c(3:length(x))])
    
    hist(val, 
         main=paste(sample," - ",gene), 
         xlab=paste("Methylation (in",length(val),"reads)"), 
         col="black",
         xlim=c(0,100),
         las=1, 
         breaks=10)
  }
} 

impHistPDF <- function (x) {
  x <- unlist(x)
  sample <- substr(x[1],19,nchar(x[1])-46)
  gene <- x[2]
  if( length(x)>2) {
    val <- as.numeric(x[c(3:length(x))])
    pdf(paste(sample,"_",gene,".pdf"),width=5,height=5,useDingbats=FALSE)
    par(pty="s")
    hist(val, 
         main=paste(sample," - ",gene), 
         xlab=paste("Methylation (in",length(val),"reads)"), 
         col="black",
         xlim=c(0,100),
         las=1, 
         breaks=10)
    dev.off()
  }
} 

filenameRAW <- file.choose() 
filenameRAW <- "/Users/meyennf/Downloads/amplicon_methylation_percentages.txt"
setwd("/Users/meyennf/Downloads/Imprint/")
data <- readLines(filenameRAW)
data <- strsplit(data, ",")
data<- strsplit(unlist (data),"\t")
x<- data[450]
for (i in 1:5)   {
  impHist(data[i])
}