#########################################################
# DESeq2 of RNASeq raw counts
#########################################################


# Get Raw Count dataset 
head(dataset)

# Make up a condition table
condition <- factor(Samples[,2])
column.data <- data.frame(source=condition)
rownames(column.data) <- Samples[,1]
column.data

#########################################################
# Load DESeq2
#########################################################
# # install.packages("DESeq2")
library("DESeq2")

# Make a DESeq data set from the counts and the design and specify which factors in the design to test

count.data.set <- DESeqDataSetFromMatrix(countData = datasetRaw, colData = column.data,design= ~ source)

# Perform the analysis.
count.data.set <- DESeq(count.data.set)


# Retrieve the full set of results.
binomial.result <- results(count.data.set,independentFiltering=FALSE)

# Remove unmeasured results
na.omit(binomial.result) -> binomial.result

# Select significant hits
binomial.result[binomial.result$padj <= 0.05,] -> significant.results

# Order the hits by p-value
significant.results[order(significant.results$padj),] -> significant.results

## Examine plot of p-values
hist(binomial.result$pvalue, breaks=50, col="grey")

## MA plot
plotMA(count.data.set, ylim=c(-1,1), cex=1)

# Write the hit names to a file
#write.table(significant.results,file="hits.txt",row.names=TRUE,col.names=NA,quote=FALSE,sep="\t")

# Get differential expression results

table(count.data.set$padj<0.05)
## Order by adjusted p-value
count.data.set <- count.data.set[order(count.data.set$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(count.data.set), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
#write.csv(resdata, file="diffexpr-results.csv")



# Plot dispersions
#png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(count.data.set, main="Dispersion plot")
#dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(count.data.set)
hist(assay(rld))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
#png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
#dev.off()


## Volcano plot with "significant" genes labeled
volcanoplot <- function (count.data.set, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(count.data.set, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(count.data.set, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(count.data.set, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(count.data.set, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(count.data.set, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
#png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2),ylim=c(0, 15))
#dev.off()







plot(all.hits$logFC,
     0-log10(all.hits$P.Value),
     xlim=c(0-max(abs(all.hits$logFC)),max(abs(all.hits$logFC))),
     pch=19,
     cex=0.4,
     main=title)

points(all.hits$logFC[all.hits$adj.P.Val < 0.05],
       (0-log10(all.hits$P.Value))[all.hits$adj.P.Val < 0.05],
       col="red"