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

#NormRow <- function(x){   (x-min(x))/(max(x)-min(x))   }

#NormDataset <- function(x) {
#  NormX <- x
#  for (i in 1:nrow(x))
#  { NormX[i,] <- cbind(x[i,c(1,2)],NormRow(2^x[i,c(-1,-2)]))  } 
#  return(NormX)
#}


DataSetName <- "RNAseq_Pipeline_All.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
datasetInputT <- datasetInput[,c(1,2,14:ncol(datasetInput))]

Transcription <- NormDataset(datasetInputT)
names(Transcription) <- c("Probe","Chromosome",rep(c("Ctrl_3CE","Ctrl_RT","PCE_3CE","PCE_RT"),each=3))

DataSetName <- "RNAseq_ActiveTranscription_All.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
datasetInputAT <- datasetInput[,c(1,2,14:ncol(datasetInput))]

ActiveTranscription <- NormDataset(datasetInputAT)
names(ActiveTranscription) <-   c("Probe","Chromosome",rep(c("Ctrl_3CE","Ctrl_RT","PCE_3CE","PCE_RT"),each=3))



par(mfrow = c(2, 1))

boxplot(Transcription[,c(-1,-2)],
        #names = c(rep("Ctrl-3CE",3),rep("Ctrl-RT",3),rep("PCE-3CE",3),rep("PCE-RT",3)),
        boxwex = 0.5, par(cex.axis=0.7),
        col = c(rep(c("light blue","orange","blue","red"),each=3)),
        main = "Exonic Transcription",
        ylab = "Normalised Expression")

boxplot(ActiveTranscription[,c(-1,-2)],
        #names = c(rep("Ctrl-3CE",3),rep("Ctrl-RT",3),rep("PCE-3CE",3),rep("PCE-RT",3)),
        boxwex = 0.5, par(cex.axis=0.7),
        col = c(rep(c("light blue","orange","blue","red"),each=3)),
        main = "Active/De novo Transcription",
        ylab = "Normalised Expression")

par(mfrow = c(2, 1))

beanplot(Transcription[,c(-1,-2)],
         #names = c(rep("Ctrl-3CE",3),rep("Ctrl-RT",3),rep("PCE-3CE",3),rep("PCE-RT",3)),
         col = "black",
         main = "Exonic Transcription",
         ylab = "Normalised Expression",
         what=c(0,1,1,0), border = NA, maxwidth=1)

beanplot(ActiveTranscription[,c(-1,-2)],
         #names = c(rep("Ctrl-3CE",3),rep("Ctrl-RT",3),rep("PCE-3CE",3),rep("PCE-RT",3)),
         col = "red",
         main = "Active/De novo Transcription",
         ylab = "Normalised Expression",
         what=c(0,1,1,0), border = NA, maxwidth=1)

AT <- cbind(ActiveTranscription[,c(1,2)],
            Ctrl_3CE = apply(ActiveTranscription[,c(3:5)],1,mean),
            Ctrl_RT = apply(ActiveTranscription[,c(6:8)],1,mean),
            PCE_3CE = apply(ActiveTranscription[,c(9:11)],1,mean),
            PCE_RT = apply(ActiveTranscription[,c(12:14)],1,mean))

TT <- cbind(Transcription[,c(1,2)],
            Ctrl_3CE = apply(Transcription[,c(3:5)],1,mean),
            Ctrl_RT = apply(Transcription[,c(6:8)],1,mean),
            PCE_3CE = apply(Transcription[,c(9:11)],1,mean),
            PCE_RT = apply(Transcription[,c(12:14)],1,mean))

par(mfrow = c(2, 1))

boxplot(TT[,c(-1,-2)],
        #names = c(rep("Ctrl-3CE",3),rep("Ctrl-RT",3),rep("PCE-3CE",3),rep("PCE-RT",3)),
        boxwex = 0.5, par(cex.axis=0.7),
        col = c(rep(c("light blue","orange","blue","red"),each=1)),
        main = "Exonic Transcription",
        ylab = "Normalised Expression")

boxplot(AT[,c(-1,-2)],
        #names = c(rep("Ctrl-3CE",3),rep("Ctrl-RT",3),rep("PCE-3CE",3),rep("PCE-RT",3)),
        boxwex = 0.5, par(cex.axis=0.7),
        col = c(rep(c("light blue","orange","blue","red"),each=1)),
        main = "Active/De novo Transcription",
        ylab = "Normalised Expression")