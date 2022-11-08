#########################################################
# Import of RNASeq
#########################################################

#Mouse DataSets Location

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/")


#########################################################
#Import combined DataSet
DataSetName <- "musPGCLC_RNASeq_OUT_20160614.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
dataset <- datasetInput[,c(1,13:ncol(datasetInput))]
remove(datasetInput)
sampleNames <- cbind(DataGroup=c(rep("Seisenberger",2),rep("Yamaguchi",4),rep("ESC",3),rep("EpiLC",3),rep("PGCLC",8)),Sample=names(dataset)[-1])

#Reduce DataSet Number
rem <- c(3:6)
dataset <- dataset[-(rem+1)]
sampleNames <- sampleNames[-rem,]
#########################################################


#########################################################
#Import all DataSets
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