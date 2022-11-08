#########################################################
# Import of Data
#########################################################

setwd("/Users/meyennf/Cambridge/Projects/Collaboration_Wolfrum/2017_PCE/2017October/")

# Define Groups and SampleNames
sampleGroups <- c("Sperm: PCE","Sperm: CTRL")
sampleNames <- make.unique(c(sampleGroups,rep(sampleGroups, each = 6)),sep = " #")[-c(1:2)]
#Samples <- data.frame(sampleGroup=rep(sampleGroups, each = 6))
#rownames(Samples) <- sampleNames
Samples <- data.frame(sampleName=sampleNames,sampleGroup=rep(sampleGroups, each = 6))
# Color4Samples <- c("#ADD8E6", "#FFA500", "#FF0000", "#00008B") # c("light blue","orange","red","dark blue")
Color2Samples <- c("brown","blue")


# Import of Promoter Methylation Data
DataSetName <- "SeqMonkExport/BS_Promoter_1_1_All.txt"
OutputName <- "Results/BS_Promoter_"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
dataset <- datasetInput[,c(14:ncol(datasetInput))][,c(7:12,1:6)] #Reorder because CTRL is first and then PCE in SeqMonkExport
#rownames(dataset) <- make.unique(datasetInput[,1])
print(colnames(dataset))
colnames(dataset) <-  sampleNames
print(colnames(dataset))


remove(datasetInput,sampleNames,sampleGroups,DataSetName)