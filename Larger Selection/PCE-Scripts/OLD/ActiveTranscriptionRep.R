#########################################################
# Import of RNASeq
#########################################################

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


NormDataset <- function(x,LogTrans = FALSE, Treshold = 0) {
  NormRow <- function(x,TresholdR = 0) {
    if (mean(t(x)) < TresholdR) return(t(rep(NA, length(x)))) 
    else return((x-min(x))/(max(x)-min(x)))
  }

  NormX <- x
  if (LogTrans) for (i in 1:nrow(x)) NormX[i,] <- cbind(x[i,c(1,2)],NormRow(2^x[i,c(-1,-2)],TresholdR = 2^Treshold))
  else for (i in 1:nrow(x)) NormX[i,] <- cbind(x[i,c(1,2)],NormRow(x[i,c(-1,-2)]))
  return(NormX)
}


DataSetName <- "RNAseq_ActiveTranscription_RepSets.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
ActiveTranscriptionRep1 <- NormDataset(datasetInput[,c(1,2,14:17)])
ActiveTranscriptionRep2 <- NormDataset(datasetInput[,c(1,2,18:21)])


DataSetName <- "RNAseq_Pipeline_RepSets.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
TranscriptionRep1 <- NormDataset(datasetInput[,c(1,2,14:17)],Treshold = 1)
TranscriptionRep2 <- NormDataset(datasetInput[,c(1,2,18:21)])


par(mfrow = c(2, 2))

boxplot(TranscriptionRep1[,c(-1,-2)],
        boxwex = 0.5, par(cex.axis=0.7),
        col = c("light blue","orange","blue","red"),
        main = "Exonic Transcription",
        ylab = "Normalised Expression")
boxplot(ActiveTranscriptionRep1[,c(-1,-2)],
        boxwex = 0.5, par(cex.axis=0.7),
        col = c("light blue","orange","blue","red"),
        main = "Active/De novo Transcription",
        ylab = "Normalised Expression")

boxplot(TranscriptionRep2[,c(-1,-2)],
        boxwex = 0.5, par(cex.axis=0.7),
        col = c("light blue","orange","blue","red"),
        main = "Exonic Transcription",
        ylab = "Normalised Expression")
boxplot(ActiveTranscriptionRep2[,c(-1,-2)],
        boxwex = 0.5, par(cex.axis=0.7),
        col = c("light blue","orange","blue","red"),
        main = "Active/De novo Transcription",
        ylab = "Normalised Expression")

par(mfrow = c(2, 2))

beanplot(TranscriptionRep1[,c(-1,-2)],
         col = "black",
         main = "Exonic Transcription",
         ylab = "Normalised Expression",
         what=c(0,1,1,0), border = NA, maxwidth=1)
beanplot(ActiveTranscriptionRep1[,c(-1,-2)],
         col = "red",
         main = "Active/De novo Transcription",
         ylab = "Normalised Expression",
         what=c(0,1,1,0), border = NA, maxwidth=1)

beanplot(TranscriptionRep2[,c(-1,-2)],
         col = "black",
         main = "Exonic Transcription",
         ylab = "Normalised Expression",
         what=c(0,1,1,0), border = NA, maxwidth=1)
beanplot(ActiveTranscriptionRep2[,c(-1,-2)],
         col = "red",
         main = "Active/De novo Transcription",
         ylab = "Normalised Expression",
         what=c(0,1,1,0), border = NA, maxwidth=1)

