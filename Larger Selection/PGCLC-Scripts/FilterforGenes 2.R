#########################################################
# Extract a number of genes
#########################################################




setwd("/Users/meyennf/Desktop/Exp/")


#GeneSets 
filenames <- list.files(pattern="*.txt", full.names=F)

species <- "hum"
dataset <- datasetHum     #Human
for (fileR in filenames) {
  Genes <- read.table(fileR,as.is=T)[,1]
  match(toupper(Genes), toupper(dataset$Probe)) -> pos # Match List, ignoring case
  data <- dataset[pos,]                                # Extract DataSubSet to new DataFrame

  write.table(t(data),paste(species,"_Selected",fileR,".csv",sep=""),
              dec = ".", sep = ",", col.names = F, row.names = T)   # Export of Data 
  }

filenames <- list.files(pattern="*.txt", full.names=F)
species <- "mus"
dataset <- datasetMus     #Mouse
for (fileR in filenames) {
  Genes <- read.table(fileR,as.is=T)[,1]
  match(toupper(Genes), toupper(dataset$Probe)) -> pos # Match List, ignoring case
  data <- dataset[pos,]                                # Extract DataSubSet to new DataFrame
  
  write.table(t(data),paste(species,"_Selected",fileR,".csv",sep=""),
              dec = ".", sep = ",", col.names = F, row.names = T)   # Export of Data 
}
