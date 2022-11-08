#########################################################
# Correlation Plot - using Seqmonk output
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/20160623/")

#DataSets 
fileName <- "CorrelationTable_musPGCLC_RNASeq.txt"
CorrData <- read.delim(fileName,header=TRUE,skip=0,fill=TRUE,dec = ".")   # Reading in data

CorrData<-CorrData[,-1]

DataGroup <- c(rep("naive mESC",3),rep("mEpiLC",3),rep("mPGCLC d4",3),
              rep("mPGCLC d6",3),rep("mPGC E9.5",2),rep("mPGC E11.5",2),
              rep("mPGC E13.5",2),rep("primed mESC",2),rep("naive mESC",2))

sampleNames <- make.unique(c(levels(as.factor(DataGroup)),DataGroup),sep=" ")[-(1:length(levels(as.factor(DataGroup))))]

#samples <- cbind(DataGroup, Sample = names(CorrData))
#samples <- cbind(DataGroup, Sample = make.unique(DataGroup,sep="-"))

colnames(CorrData) <- sampleNames
rownames(CorrData) <- sampleNames


# CorrMatrix  <- data.matrix(CorrData)
# 
# tCorrMatrix <- t(data.matrix(CorrData)[ncol(CorrData):1,])
# image(tCorrMatrix,
#       main="mouse RNASeq Data",
#       
#       xlab="x column", ylab="y column",
#       col=colorRampPalette(c("red", "yellow"))(n = 100))
# 
# 
# #install.packages("corrplot")
# library(corrplot)
# corrplot(CorrMatrix,
#          type="upper",tl.col="black",
#          method="color",col=colorRampPalette(c("red", "yellow","blue","green"))(n = 100))
# 
# 
# #install.packages("corrgram")  
# library(corrgram)
# corrgram(CorrData, order=F,
#          main="mouse RNASeq Data")





library(reshape2)
melted_CorrMatrix <- melt(data.matrix(CorrData))    # Convert to Format for ggplot
melted_CorrMatrix$Var1 <- factor(melted_CorrMatrix$Var1, levels=sampleNames) # flipping the order of the y-axis
melted_CorrMatrix$Var2 <- factor(melted_CorrMatrix$Var2, levels=rev(sampleNames)) # flipping the order of the y-axis
head(melted_CorrMatrix)


library(ggplot2)

#pdf(paste(fileName,"pdf",sep="."),width=10,height=10,useDingbats=FALSE)
par(cex.lab=1, cex.axis=0.7, cex.sub=0.5)

ggheatmap <- ggplot(data = melted_CorrMatrix, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  ggtitle(fileName) +
  scale_fill_gradientn(colours=c("blue","white","red"),
                       values=c(0,0.5,1),
                       na.value = "white",
                       name="Correlation") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1)) +
  coord_fixed()+
  coord_equal() 

print(ggheatmap)

#dev.off()


