
# Identify the MspI recongnition sites for each chromosomal entry


### Load the needed libraries
# Biostrings is needed for pattern identification
# BSgenome.Mmusculus.UCSC.mm10 is the mouse genome
# plyr and reshape2 are needed for manipulating the data format

#biocLite("Biostrings")
#install.packages("scales")

library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)


# Only consider standard chromosomes
chromNames <- names(Mmusculus)[!grepl("_", names(Mmusculus))]

# Identify all CG sites in the genome
print("Finding all CpG positions")
CGs <- lapply(chromNames, function(x) start(matchPattern("CG", Mmusculus[[x]])))
names(CGs) <- chromNames
#str(CGs)

# Generate a dataframe with the length of MspI digested fragments and count CGs in each fragment
mdf=data.frame()
for (i in 1:length(chromNames)){
  print(paste("Processing ",seqnames(Mmusculus)[i], sep=""))
  m <- matchPattern("CCGG", Mmusculus[[i]])
  
  # Start and End position of each MspI site
  starts <- start(gaps(m))-4
  starts <- replace(starts, starts == -3, 0)
  ends <- end(gaps(m))
  
  # Count the number of CGs in each MspI fragment
  CG_onChrom <- unlist(CGs[chromNames[i]])
  CpGNum <- hist(CG_onChrom, c(0,ends), plot = FALSE)$counts
  
  temp_df<-data.frame(chr=seqnames(Mmusculus)[i],start=starts,end=ends,CpGs=CpGNum) #actually end = ends
  mdf<-rbind(mdf,temp_df)
}

remove(m,temp_df,i,starts,ends,CG_onChrom,CpGNum)

# Extract the digested fragment length
mdf$width=mdf$end-mdf$start

# Keep only the fragments in the specified range

ml<-mdf[mdf$width>=40&mdf$width<600,]

countFragments <- ddply(ml,.(width), nrow)
countCpGs <- ddply(ml,.(width),summarise, sum = sum(CpGs), mean = mean(CpGs),max =max(CpGs),median=median(CpGs))
  
countCpGs$meanCorr <- (countCpGs$mean ) / countCpGs$width

# Create a plot of the frequency of the fragment lengths (y-axis is logarithmic)
p<-ggplot(countFragments,aes(x=width, y=V1))+geom_line()
p+scale_y_continuous(trans=log2_trans())

p<-ggplot(countCpGs,aes(x=width, y=sum))+geom_line()
p+scale_y_continuous(trans=log2_trans())

p<-ggplot(countCpGs,aes(x=width, y=mean))+geom_line()
p+scale_y_continuous(trans=log2_trans())

p<-ggplot(countCpGs,aes(x=width, y=meanCorr))+geom_line()
p+scale_y_continuous(trans=log2_trans())

p<-ggplot(countCpGs,aes(x=width, y=max))+geom_line()
p+scale_y_continuous(trans=log2_trans())

p<-ggplot(countCpGs,aes(x=width, y=median))+geom_line()
p+scale_y_continuous(trans=log2_trans())
