########

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/Composition/")

KMER <- c("AA","AC","AG","AT","CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

# Import Human Genome Files
genomeHumFiles <- list.files(path = "Genome", pattern="Homo_sapiens*", full.names=F)
genomeHum <- rbind (EMPTY=numeric(16)); colnames(genomeHum) <- KMER
for (filenameRAW in genomeHumFiles) {
  #  filenameRAW <- "Homo_sapiens.GRCh37.55.dna.chromosome.1_compter.txt"
  Input <- read.delim(paste("Genome/",filenameRAW,sep=""))[,2]
  genomeHum <- rbind(genomeHum,Input)
  rownames(genomeHum)[nrow(genomeHum)] <- paste0("Chr",unlist(strsplit(unlist(strsplit(filenameRAW,".", fixed = TRUE))[6],"_compter")))
} 
genomeHum <- genomeHum[-1,]
hist(genomeHum[,7]); HumCG <- apply(t(genomeHum[,7]),1,mean) # Calculate Human CG enrichment average

# Import Mouse Genome Files
genomeMusFiles <- list.files(path = "Genome", pattern="Mus_musculus.*", full.names=F)
genomeMus <-  rbind (EMPTY=numeric(16)); colnames(genomeMus) <- KMER
for (filenameRAW in genomeMusFiles) {
  Input <- read.delim(paste("Genome/",filenameRAW,sep=""))[,2]
  genomeMus <- rbind(genomeMus,Input)
  rownames(genomeMus)[nrow(genomeMus)] <- paste0("Chr",unlist(strsplit(unlist(strsplit(filenameRAW,".", fixed = TRUE))[6],"_compter")))
} 
genomeMus <- genomeMus[-1,]
hist(genomeMus[,7]); MusCG <- apply(t(genomeMus[,7]),1,mean) # Calculate Mouse CG enrichment average

remove(Input,genomeMusFiles,genomeHumFiles,filenameRAW) # Cleanup




require(grid)
require(scales)


classes <- c(rep("NULL",7),"numeric",rep("NULL",9)) #Read only CG
CG <- read.table(paste("HUM/",filenameRAW,sep=""),header=TRUE,skip=0,fill=TRUE,dec = ".",colClasses = classes) - HumCG 
beanplot(CG, what=c(1,1,1,0),
         main = paste("CG enrichment\nnormalized to the genome average\n",filenameRAW),wd=1,bw="nrd0",
         col = "red")


HumFiles <- list.files(path = "HUM", pattern="*txt", full.names=F)
#classes <- c("character", rep("NULL",6),"numeric",rep("NULL",9)) #Read Position and CG
classes <- c(rep("NULL",7),"numeric",rep("NULL",9)) #Read only CG
require(beanplot) #install.packages("beanplot")
pdf("HUMAN_BeanPlots.pdf")
for (filenameRAW in HumFiles) {
  CG <- read.table(paste("HUM/",filenameRAW,sep=""),header=TRUE,skip=0,fill=TRUE,dec = ".",colClasses = classes) - HumCG 
  beanplot(CG, what=c(1,1,1,0),
           main = paste("CG enrichment\nnormalized to the genome average\n",filenameRAW),wd=1,bw="nrd0",
           col = "red")
} 
dev.off()

filenameRAW <- "/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/Composition/MUS/T_musBSData_compter.txt"

MusFiles <- list.files(path = "MUS", pattern="*txt", full.names=F)
#classes <- c("character", rep("NULL",6),"numeric",rep("NULL",9)) #Read Position and CG
classes <- c(rep("NULL",7),"numeric",rep("NULL",9)) #Read only CG
require(beanplot) #install.packages("beanplot")
pdf("Mouse_BeanPlots.pdf")
for (filenameRAW in MusFiles) {
  CG <- read.table(paste("MUS/",filenameRAW,sep=""),header=TRUE,skip=0,fill=TRUE,dec = ".",colClasses = classes) - MusCG 
  CG <- read.table(filenameRAW,header=TRUE,skip=0,fill=TRUE,dec = ".",colClasses = classes) - MusCG 
  
  beanplot(CG, what=c(1,1,1,0),
           main = paste("CG enrichment\nnormalized to the genome average\n",filenameRAW),wd=1,bw="nrd0",
           col = "blue")
} 
dev.off()