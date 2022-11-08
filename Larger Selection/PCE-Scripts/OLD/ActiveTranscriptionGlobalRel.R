setwd("/Users/meyennf/Cambridge/Projects/Collaborations/Christian Wolfrum/PCE/RNAseq/")

NormDataset <- function(x,LogTrans = FALSE, Treshold = 0) {
  
  NormRow <- function(x,TresholdR = 0) {
    if (mean(t(x)) < TresholdR) return(t(rep(NA, length(x)))) 
    else return((x-mean(t(x)))/sd(t(x)))
    # else return((x-min(x))/(max(x)-min(x)))
  }
  
  NormX <- x
  if (LogTrans) for (i in 1:nrow(x)) NormX[i,] <- cbind(x[i,c(1,2)],NormRow(2^x[i,c(-1,-2)],TresholdR = 2^Treshold))
  else for (i in 1:nrow(x)) NormX[i,] <- cbind(x[i,c(1,2)],NormRow(x[i,c(-1,-2)]))
  return(NormX)
}

#NormRow <- function(x){   (x-min(x))/(max(x)-min(x))   }

#NormDataset <- function(x) {
#  NormX <- x
#  for (i in 1:nrow(x))
#  { NormX[i,] <- cbind(x[i,c(1,2)],NormRow(2^x[i,c(-1,-2)]))  } 
#  return(NormX)
#}


datasetInput <- read.csv("RNAseq_Pipeline_All.txt",sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
datasetInputT <- datasetInput[,c(1,2,14:ncol(datasetInput))]

datasetInput <- read.csv("RNAseq_ActiveTranscription_All.txt",sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
datasetInputAT <- datasetInput[,c(1,2,14:ncol(datasetInput))]


Transcription <- NormDataset(datasetInputT)
names(Transcription) <- c("Probe","Chromosome",rep(c("Ctrl_3CE","Ctrl_RT","PCE_3CE","PCE_RT"),each=3))


ActiveTranscription <- NormDataset(datasetInputAT)
names(ActiveTranscription) <-   c("Probe","Chromosome",rep(c("Ctrl_3CE","Ctrl_RT","PCE_3CE","PCE_RT"),each=3))
