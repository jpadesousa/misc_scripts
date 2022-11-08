#########################################################
# Import of RNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/")

DataSetName <- "musPGCLC_RNASeq_OUT.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

dataset <- datasetInput[,c(1,13:ncol(datasetInput))]

sampleNames <- cbind(c(rep("Merged",4),rep("Unique",16)),names(dataset)[-1])

