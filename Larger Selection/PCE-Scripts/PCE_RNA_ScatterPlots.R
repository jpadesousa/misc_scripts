
# Import of RNASeq - SignificantGeneList (exported from SeqMonk)
source('~/Cambridge/Projects/Collaboration_Wolfrum/2017_PCE/2017October/Scripts/PCE_RNA_PreProcessing.R')

GeneListFolderMax <- "Results_RNA/7_RNA_Lists/1-MaxStringency/"
GeneListFileNamesMax <- list.files(GeneListFolderMax, pattern=".txt")
GeneListFolderMed <- "Results_RNA/7_RNA_Lists/2-MedStringency/"
GeneListFileNamesMed <- list.files(GeneListFolderMed, pattern=".txt")
ListNames <- unlist(strsplit(GeneListFileNamesMed,"_Med.txt"))

# Gene Lists to Highlight
GeneListCol <- list()
GeneListLabel <- list()
for (i in 1:4) {
  print(ListNames[i])
  GeneListCol[[ListNames[i]]] <- read.csv(paste(GeneListFolderMed,GeneListFileNamesMed[i],sep=""),
                                          sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
  GeneListLabel[[ListNames[i]]] <- read.csv(paste(GeneListFolderMax,GeneListFileNamesMax[i],sep=""),
                                            sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
}

# Data to plot
dataPlotList <- list()
# 1 Ctrl: 3CE vs RT
dataPlotList[[ListNames[1]]] <- cbind(Probe = rownames(dataset), 'BAT: Ctrl-CE' = apply(dataset[,c(7:12)],1,mean),
                                      'BAT: Ctrl-RT' =  apply(dataset[,c(19:24)],1,mean))
# 2 PCE: 3CE vs RT
dataPlotList[[ListNames[2]]] <- cbind(Probe = rownames(dataset), 'BAT: PCE-CE' =  apply(dataset[,c(1:6)],1,mean),
                                      'BAT: PCE-RT' =  apply(dataset[,c(13:18)],1,mean))
# 3 PCE vs CTRL (at 3CE)
dataPlotList[[ListNames[3]]] <- cbind(Probe = rownames(dataset), 'BAT: PCE-CE' =  apply(dataset[,c(1:6)],1,mean),
                                      'BAT: Ctrl-CE' =  apply(dataset[,c(7:12)],1,mean))
# 4 PCE vs CTRL (at RT)
dataPlotList[[ListNames[4]]] <-  cbind(Probe = rownames(dataset), 'BAT: PCE-RT' =  apply(dataset[,c(13:18)],1,mean),
                                       'BAT: Ctrl-RT' =  apply(dataset[,c(19:24)],1,mean))



for (i in 1:4) {

x <- as.numeric (dataPlotList[[i]][,2])
y <- as.numeric (dataPlotList[[i]][,3])

KeyGenes <- c("Ucp1","Adrb3","Bmp8b", "Mfsd2a", "Fndc5", "Alpl", "Ckm", "Ckmt2", "Gpr64", "Gyk","Cidea","Dio2","Pparg","Fabp4")
posListKey <- (match(dataPlotList[[i]][,1],setdiff(KeyGenes,GeneListCol[[i]][,1])) > 0)
#cbind(KeyGenes,match(KeyGenes,dataPlotList[[i]][,1]))

posListCol <- (match(dataPlotList[[i]][,1],GeneListCol[[i]][,1]) > 0)
posListLabels <- (match(dataPlotList[[i]][,1],
                        union(GeneListLabel[[i]][,1],intersect(KeyGenes,GeneListCol[[i]][,1]))) > 0)

df <- data.frame(x,y)

RSQ <- round(cor.test(x, y)[4][[1]],3)
      
## Use densCols() output to get density at each point
      
x <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L
      
## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                                  "#FCFF00", "#FF9400", "#FF3100"))(256)

df$col <- cols[df$dens]
xmin <- floor(as.numeric(min(x,y)))
xmax <- ceiling(as.numeric(max(x,y)))
 

## Create plot for R 
# par(pty="s")
# plot(y~x, data=df[order(df$dens),], # col=cols, # Reordering rows so that densest points are plotted on top
#      pch=19, cex=.2,type ="p", 
#      xlim=c(xmin,xmax),ylim=c(xmin,xmax),
#      axes = T,frame.plot = F,asp=1,
#      ylab = colnames(dataPlotList[[i]])[3] ,
#      xlab = colnames(dataPlotList[[i]])[2])
# 
# points(y~x, data=df[posListCol,], pch=19, col="red", cex=.2,type ="p")
# 
# text(y~x, data=df[posListLabels,], labels=dataPlotList[[i]][posListLabels,1], cex= 0.6,pos=2*(1+as.integer(df$x / df$y > 1)))
# text(xmax-4,xmin,paste("R =",round(RSQ,3)),cex= 0.9)



## Create plot and capture as image and use to reduce reendering in the PDF
require(grid)
dev.new(noRStudioGD=T,height=7,width=7)
par(mar=c(0,0,0,0),mai=c(0,0,0,0),pty="s") 
plot(y~x, data=df[order(df$dens),],col="dark grey", # col=cols, # Reordering rows so that densest points are plotted on top
     pch=19, cex=.2,type ="p", 
     xlim=c(xmin,xmax),ylim=c(xmin,xmax),
     axes = F,frame.plot = F,asp=1,#bty ="n",
     xaxs = "i", yaxs = "i",xlab = "",ylab="")
points(y~x, data=df[posListKey,], pch=19, col="blue", cex=.2,type ="p")
points(y~x, data=df[posListCol,], pch=19, col="red", cex=.2,type ="p")
capture <- grid.cap()
dev.off()
    
#Create PDF with plot + Frame
PDFName <- paste("Results_RNA/RNA_Scatter_",ListNames[i],".pdf",sep="")
pdf(PDFName,width=10,height=10,useDingbats=FALSE)
par(pty="s")
plot(1,type="n",
     xlim=c(xmin,xmax),ylim=c(xmin,xmax),
     axes = T,frame.plot = F,asp=1,
     ylab = paste(colnames(dataPlotList[[i]])[3],"(log2 RPM)") ,
     xlab = paste(colnames(dataPlotList[[i]])[2],"(log2 RPM)") )

rasterImage(capture,xmin,xmin,xmax,xmax)
text(y~x, data=df[posListKey,], labels=dataPlotList[[i]][posListKey,1], cex= 0.6,pos=2*(1+as.integer(df$x / df$y > 1)),col="blue")
text(y~x, data=df[posListLabels,], labels=dataPlotList[[i]][posListLabels,1], cex= 0.6,pos=2*(1+as.integer(df$x / df$y > 1)))

text(xmax,xmin,paste("R =",round(RSQ,3)),cex= 0.9,pos=2)
legend(xmin,xmax, legend=c("Sign. Differentially Expressed", "Selected Genes of Interest (not DE)"),
       col=c("red", "blue"), lty=1,lwd=2, cex=0.6,bty="n")


dev.off()
}

