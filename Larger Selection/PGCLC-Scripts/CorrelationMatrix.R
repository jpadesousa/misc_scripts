#########################################################
# Correlation Plot - using Seqmonk output
#########################################################

setwd("/Users/meyennf/Documents/Cambridge/Projects/PGCLC/2016Data/musRNASeq/")

#DataSets 
fileName <- "Correlation_mPGC_RNA.txt"
CorrData <- read.delim(fileName,header=TRUE,skip=0,fill=TRUE,dec = ".")   # Reading in data

#rownames(CorrData) <- CorrData[,1]
CorrData<-CorrData[,-1]

c("E14_ser_1","E14_ser_2","E14_ser_3",
  "E14_2i_1","E14_2i_2","E14_2i_3",
  "E14_EpiLCd2_2","E14_EpiLCd2_3","E14_EpiLCd2_1",
  "mPGCLC d4 E2","mPGCLC d4 F2","mPGCLC d4 G2","mPGCLC d4 H2",
  "mPGCLC d6 A5","mPGCLC d6 B5","mPGCLC d6 C5","mPGCLC d6 D5") -> SAMPLENAMES


colnames(CorrData) <- SAMPLENAMES
rownames(CorrData) <- SAMPLENAMES

CorrMatrix <- scale(CorrMatrix)



CorrMatrix  <- data.matrix(CorrData)
tCorrMatrix <- t(data.matrix(CorrData)[ncol(CorrData):1,])

image(tCorrMatrix,
      main="mouse RNASeq Data",
      
      xlab="x column", ylab="y column",
      col=colorRampPalette(c("red", "yellow"))(n = 100))


#install.packages("corrplot")
library(corrplot)
corrplot(CorrMatrix,
         type="upper",tl.col="black",
         method="color",col=colorRampPalette(c("red", "yellow","blue","green"))(n = 100))


#install.packages("corrgram")  
library(corrgram)
corrgram(CorrData, order=F,
         main="mouse RNASeq Data")





library(reshape2)
melted_CorrMatrix <- melt(data.matrix(CorrData))    # Convert to Format for ggplot
melted_CorrMatrix$Var1 <- factor(melted_CorrMatrix$Var1, levels=SAMPLENAMES) # flipping the order of the y-axis
melted_CorrMatrix$Var2 <- factor(melted_CorrMatrix$Var2, levels=rev(SAMPLENAMES)) # flipping the order of the y-axis
head(melted_CorrMatrix)


library(ggplot2)
ggheatmap <- ggplot(data = melted_CorrMatrix, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradientn(colours=c("blue","white","red"),
                       values=c(0,0.5,1),
                       #  guide = "none",
                       na.value = "white",
                       name="Correlation") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1)) +
  coord_fixed()+
  coord_equal() 

print(ggheatmap)
