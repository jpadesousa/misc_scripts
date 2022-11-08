#########################################################
# Plotting a heat map
#########################################################

#install.packages("gplots")
#install.packages("RColorBrewer")

library(gplots)
library(RColorBrewer)

#########################################################
# Reading in data and transform it into matrix format
#########################################################

mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- data[,1]                # assign row names

#########################################################
# Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_Col_palette <- colorRampPalette(c("green", "red"))(n = 299)

ColorSet <- rainbow(length(unique(sampleNames[,1])))
names(ColorSet) <- unique(sampleNames[,1])
ColorLabels=ColorSet[sampleNames[,1]]

# Clustering

hclustfunc <- function(x) hclust(x, method="complete")
#hclustfunc <- function(x) hclust(x, method="ward.D")

distfunc <- function(x) dist(x,method="euclidean")
#distfunc <- function(x) as.dist((1-cor(t(x)))/2) # Pearson
#pdf(paste(DataSetName,"pdf",sep="."),width=10,height=10,useDingbats=FALSE)
#par(cex.lab=1, cex.axis=0.7, cex.sub=0.5)
heatmap.2(mat_data,
          main = paste("Correlation in",DataSetName), # heat map title
          density="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(8,8),     # widens margins around plot
          col = my_Col_palette,       # use on color palette defined earlier
          dendrogram = "col",     # only draw a row or col dendrogram
          scale = "row",
          #RowSideColors = rep("red",nrow(mat_data)),
          ColSideColors = ColorLabels,
          hclustfun = hclustfunc,distfun = distfunc,
          srtCol=45,cexCol=0.7,offsetCol=-0.5)            

legend("bottomleft", legend = unique(sampleNames[,1]), 
       col = ColorSet, lty= 1, lwd = 5, cex=0.7)          
#dev.off()

