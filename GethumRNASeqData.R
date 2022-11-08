#########################################################
# Import of RNASeq
#########################################################

#Human DataSets Location
setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/humRNASeq/")

#########################################################
fileName <- "20160624/Expression_humPGCLC_RNASeq.txt"
dataset <- read.csv(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
dataset <- dataset[,c(1,13:ncol(dataset))]

#########################################################
#Reduce DataSet Number
rem <- c(7:10,22:24)
dataset <- dataset[-(rem+1)]
sampleNames <- sampleNames[-rem,]
#########################################################

DataGroup <- c(rep("naive hESC",3),rep("primed hESC",3),rep("hPGCLC",5),rep("hPGC",6),rep("hSoma",3))

#sampleNames <- make.unique(c(levels(as.factor(DataGroup)),DataGroup),sep=" ")[-(1:length

sampleNames <- c("naive hESC 1","naive hESC 2","naive hESC 3",
                 "primed hESC 1","primed hESC 2","primed hESC 3",
                 "hPGCLC d4 1","hPGCLC d4 2","hPGCLC d5",
                 "hPGCLC d8","hPGCLC d12","hPGC Wk5.5",
                 "hPGC Wk7F 1","hPGC Wk7F 2","hPGC Wk7F 3",
                 "hPGC Wk9F 1","hPGC Wk9F 2",
                 "hSoma 1","hSoma 2","hSoma 3")


#Samples <- cbind(DataGroup,names(dataset)[-1])
Samples <- cbind(DataGroup,sampleNames)

#sampleNames <- cbind(DataGroup=c(rep("hPrimed",3),rep("hNaive",3),rep("hEpiLC",4),rep("hPGCLC",5),rep("hPGC",6),rep("hSoma",3)),Sample=names(dataset)[-1])



#########################################################
#Import all Human DataSets
filenames <- list.files(path = "Data", pattern="*.txt", full.names=F)
for (input in filenames) {
  datasetInput <- read.csv(paste("Data/",input,sep=""),sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
  if (input == filenames[1])
    dataset <- datasetInput[,c(1,13:ncol(datasetInput))]
  else {
    match(dataset[,1], datasetInput[,1]) -> pos
    if (all(dataset[,1] == datasetInput[pos,1]))
      dataset <- cbind(dataset,datasetInput[pos,c(13:ncol(datasetInput))])
  }
}
remove(input,pos,datasetInput)
sampleSet <- read.csv("SampleNames.txt",sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
names(dataset) <-   c("Probe",sampleSet[,1])
sampleNames <- sampleSet
########################################################