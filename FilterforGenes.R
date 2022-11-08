#########################################################
# Extract a number of genes
#########################################################
#setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/")
Genes <- read.table("GenesK9.txt",as.is=T)[,1]


#dataset <- datasetHum
#dataset <- datasetMus

# Match List, ignoring case
match(toupper(Genes), toupper(RNAseqdataset$Probe)) -> pos




#cbind(Genes,pos,dataset[pos,1])
#match("DKFZp434M202", data[,1])

#Extract DataSubSet to new DataFrame
data <- dataset[pos,]
remove(pos)


# Remove Rows with NAs
# data <- na.omit(data)

#########################################################
# Export of Data 
#########################################################
setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/Figures/v3_July/Figure4/PIWI_Expression/")
#output <- rbind(c("DataSet",sampleNames[,2]),data)
write.table(t(data),
            "PiwiGenes_Human.txt",
#            paste("SelectedGenes_",DataSetName,".txt",sep=""),
            sep="\t",dec = ".",row.names = T,col.names = F)
#remove(output)
#########################################################
# Extract samples from specific DataSets
#########################################################


# FocusSet <- c("Gkountela","Irie","hPGCLC")
# 
# !is.na(match(sampleNames[,2],FocusSet)) -> pos
# 
# #Extract DataSubSet to new DataFrame & Adjust sampleNames
# data <- dataset[,pos]
# sampleNames <- sampleSet[pos,]
# 
# remove(pos)
