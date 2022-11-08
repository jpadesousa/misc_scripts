#########################################################
# Import of RNASeq
#########################################################

setwd("/Users/meyennf/Cambridge/Projects/Collaboration_Wolfrum/2017_PCE/2017October/")

# Define Groups and SampleNames
sampleGroups <- c("BAT: PCE-CE","BAT: CTRL-CE","BAT: PCE-RT","BAT: CTRL-RT")
sampleNames <- make.unique(c(sampleGroups,rep(sampleGroups, each = 6)),sep = " #")[-c(1:4)]
#Samples <- data.frame(sampleGroup=rep(sampleGroups, each = 6))
#rownames(Samples) <- sampleNames
Samples <- data.frame(sampleName=sampleNames,sampleGroup=rep(sampleGroups, each = 6))

Color4Samples <- c("#ADD8E6", "#FFA500", "#FF0000", "#00008B") # c("light blue","orange","red","dark blue")


# Import of RNASeq - Norm Counts
DataSetName <- "SeqMonkExport/RNAseq_V2_NormCounts.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
dataset <- datasetInput[,c(14:ncol(datasetInput))]
rownames(dataset) <- make.unique(datasetInput[,1])
colnames(dataset) <-  sampleNames

# Import of RNASeq - Raw Counts
DataSetName <- "SeqMonkExport/RNAseq_V2_RawCounts.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
datasetRaw <- datasetInput[,c(14:ncol(datasetInput))]
rownames(datasetRaw) <- make.unique(datasetInput[,1])
colnames(datasetRaw) <- sampleNames



remove(datasetInput,sampleNames,sampleGroups,DataSetName)

