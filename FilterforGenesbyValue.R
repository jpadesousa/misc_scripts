#########################################################
# Filter for Geneexpressison: Comparison is d0 vs d10!
#########################################################



filteredRNAseqdata <- RNAseqdataset [1:100,]

# Remove genes with no information (remove NA rows)


# 1) Remove genes with p-value > 0.00005

filteredRNAseqdata <- filteredRNAseqdata[filteredRNAseqdata$P.value < 0.00005,]

# 2) Remove genes with very low expression in d0 and d10 (expression in both <= 25 reads remove)

# filteredRNAseqdata <- filteredRNAseqdata[!(filteredRNAseqdata$D0mean < 25 & filteredRNAseqdata$D10mean < 25),]


# List of HIGH d0 - HIGH d10 Genes (high: > 500 reads )
HighHigh <- filteredRNAseqdata[(filteredRNAseqdata$D0mean > 500 & filteredRNAseqdata$D10mean > 500),]

# List of HIGH d0 - Low d10 Genes (high: > 500 reads )
HighLow <- filteredRNAseqdata[(filteredRNAseqdata$D0mean > 500 & filteredRNAseqdata$D10mean < 100),]

# List of Low d0 - Low d10 Genes (high: > 500 reads )
LowLow <- filteredRNAseqdata[(filteredRNAseqdata$D0mean < 100 & filteredRNAseqdata$D10mean < 100),]

# List of Low d0 - High d10 Genes (high: > 500 reads )
LowHigh <- filteredRNAseqdata[(filteredRNAseqdata$D0mean < 100 & filteredRNAseqdata$D10mean > 500),]



Genes <- LowLow$Probe
