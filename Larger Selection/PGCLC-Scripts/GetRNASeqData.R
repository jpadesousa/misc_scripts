#########################################################
# Import of huRNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/humRNASeq/20160624_2/")

fileName <- "Expression_humPGCLC_RNASeq.txt"
dataset <- read.csv(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

datasetHum <- dataset[,c(1,13:ncol(dataset))]

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
SamplesHum <- cbind(DataGroup,sampleNames)
names(datasetHum) <- c("Probe",sampleNames)
remove(DataGroup,fileName,sampleNames,dataset)
#########################################################
#Remove Samples
SamplesHum <- SamplesHum[-c(12:17,18:20 ),]
datasetHum <- datasetHum[,-(1+c(12:17,18:20))]


#########################################################
# Import of mus RNASeq
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/20160624/")

fileName <- "Expression_musPGCLC_RNASeq.txt"
dataset <- read.csv(fileName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)

datasetMus <- dataset[,c(1,13:ncol(dataset))]

DataGroup <- c(rep("naive mESC",3),rep("mEpiLC",3),rep("mPGCLC d4",3),
               rep("mPGCLC d6",3),rep("mPGC E9.5",2),rep("mPGC E11.5",2),
               rep("mPGC E13.5",2),rep("primed mESC",2),rep("naive mESC",2))

sampleNames <- make.unique(c(levels(as.factor(DataGroup)),DataGroup),sep=" ")[-(1:length(levels(as.factor(DataGroup))))]

#Samples <- cbind(DataGroup,names(dataset)[-1])

SamplesMus <- cbind(DataGroup,sampleNames)
names(datasetMus) <- c("Probe",sampleNames)
remove(DataGroup,fileName,sampleNames,dataset)
#########################################################

#Remove Samples
SamplesMus <- SamplesMus[-c(13:18,19:22),]
datasetMus <- datasetMus[,-(1+c(13:18,19:22))]
#########################################################

