# #########################################################  
# GO Analysis using gProfileR
# #########################################################  

# Get Gene Lists
source('~/Cambridge/Projects/Collaboration_Wolfrum/2017_PCE/2017October/Scripts/PCE_RNA_PreProcessing.R')
setwd("/Users/meyennf/Cambridge/Projects/Collaboration_Wolfrum/2017_PCE/2017October/")

DataSetName <- "SeqMonkExport/RNABS_GenesOverlappingHyperPCE.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
DMRGeneListHyper <- datasetInput[,1]

DataSetName <- "SeqMonkExport/RNABS_GenesOverlappingHypoPCE.txt"
datasetInput <- read.csv(DataSetName,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",as.is=T)
DMRGeneListHypo <- datasetInput[,1]

GeneListFolderMed <- "Results_RNA/7_RNA_Lists/2-MedStringency/"
GeneListFileNamesMed <- list.files(GeneListFolderMed, pattern=".txt")
ListNames <- unlist(strsplit(GeneListFileNamesMed,"_Med.txt"))
remove(datasetInput,DataSetName,Samples,datasetRaw)

GeneList <- list()
ToGOListALL <- list()
ToGOListUP <- list()
ToGOListDOWN <- list()

for (i in 1:4) {
  GeneList[[ListNames[i]]] <- read.csv(paste(GeneListFolderMed,GeneListFileNamesMed[i],sep=""),
                                       sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
  ToGOListALL[[ListNames[i]]] <- as.character(GeneList[[i]][,1])
}

# 1 Ctrl: 3CE vs RT: CTRL_3CEvsRT
ToGOListUP$CTRL_3CEvsRT   <-  as.character(GeneList$CTRL_3CEvsRT[(GeneList$CTRL_3CEvsRT$CTRL_3CE - 
                                                                    GeneList$CTRL_3CEvsRT$CTRL_RT > 0),1])
ToGOListDOWN$CTRL_3CEvsRT <-  as.character(GeneList$CTRL_3CEvsRT[(GeneList$CTRL_3CEvsRT$CTRL_3CE - 
                                                                    GeneList$CTRL_3CEvsRT$CTRL_RT < 0),1])

# 2 PCE: 3CE vs RT: PCE_3CEvsRT
ToGOListUP$PCE_3CEvsRT   <-  as.character(GeneList$PCE_3CEvsRT[(GeneList$PCE_3CEvsRT$PCE_3CE - 
                                                                  GeneList$PCE_3CEvsRT$PCE_RT > 0),1])
ToGOListDOWN$PCE_3CEvsRT <-  as.character(GeneList$PCE_3CEvsRT[(GeneList$PCE_3CEvsRT$PCE_3CE - 
                                                                  GeneList$PCE_3CEvsRT$PCE_RT < 0),1])

# 3 PCE vs CTRL (at 3CE): PCEvsCTRL_3CE
ToGOListUP$PCEvsCTRL_3CE   <-  as.character(GeneList$PCEvsCTRL_3CE[(GeneList$PCEvsCTRL_3CE$PCE_3CE - 
                                                                      GeneList$PCEvsCTRL_3CE$CTRL_3CE > 0),1])
ToGOListDOWN$PCEvsCTRL_3CE <-  as.character(GeneList$PCEvsCTRL_3CE[(GeneList$PCEvsCTRL_3CE$PCE_3CE - 
                                                                      GeneList$PCEvsCTRL_3CE$CTRL_3CE < 0),1])

# 4 PCE vs CTRL (at RT): PCEvsCTRL_RT
ToGOListUP$PCEvsCTRL_RT   <-  as.character(GeneList$PCEvsCTRL_RT[(GeneList$PCEvsCTRL_RT$PCE_RT - 
                                                                    GeneList$PCEvsCTRL_RT$CTRL_RT > 0),1])
ToGOListDOWN$PCEvsCTRL_RT <-  as.character(GeneList$PCEvsCTRL_RT[(GeneList$PCEvsCTRL_RT$PCE_RT - 
                                                                    GeneList$PCEvsCTRL_RT$CTRL_RT < 0),1])



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

i=1

for (i in 1:4) {
  x <- as.numeric (dataPlotList[[i]][,2])
  y <- as.numeric (dataPlotList[[i]][,3])
 
  ALL_Hyper <- intersect(ToGOListALL[[i]],DMRGeneListHyper)
#  UP_Hyper <-  intersect(ToGOListUP[[i]],DMRGeneListHyper)
#  DOWN_Hyper <- intersect(ToGOListDOWN[[i]],DMRGeneListHyper)
  
  ALL_Hypo <- intersect(ToGOListALL[[i]],DMRGeneListHypo)
#  UP_Hypo <- intersect(ToGOListUP[[i]],DMRGeneListHypo)
#  DOWN_Hypo <- intersect(ToGOListDOWN[[i]],DMRGeneListHypo)
  
  posList_ALL_Hyper <- (match(dataPlotList[[i]][,1],ALL_Hyper) > 0)
#  posList_UP_Hyper <- (match(dataPlotList[[i]][,1],UP_Hyper) > 0)
#  posList_DOWN_Hyper <- (match(dataPlotList[[i]][,1],DOWN_Hyper) > 0)
  
  posList_ALL_Hypo <- (match(dataPlotList[[i]][,1],ALL_Hypo) > 0)
#  posList_UP_Hypo <- (match(dataPlotList[[i]][,1],UP_Hypo) > 0)
#  posList_DOWN_Hypo <- (match(dataPlotList[[i]][,1],DOWN_Hypo) > 0)
 
  posListLabels <- posList_ALL_Hyper | posList_ALL_Hypo 
  
#  KeyGenes <- c("Ucp1","Adrb3","Bmp8b", "Mfsd2a", "Fndc5", "Alpl", "Ckm", "Ckmt2", "Gpr64", "Gyk","Cidea","Dio2","Pparg","Fabp4")
#  cbind(KeyGenes,match(KeyGenes,dataPlotList[[i]][,1]))
  
#  posListKey <- (match(dataPlotList[[i]][,1],KeyGenes) > 0)
#  posListLabels <- posListLabels | posListKey
  
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
  #     pch=19, cex=.2,type ="p",
  #     xlim=c(xmin,xmax),ylim=c(xmin,xmax),
  #     axes = T,frame.plot = F,asp=1,
  #     ylab = colnames(dataPlotList[[i]])[3] ,
  #     xlab = colnames(dataPlotList[[i]])[2])
  # 
  # points(y~x, data=df[posList_ALL_Hyper,], pch=19, col="red", cex=.2,type ="p")
  # points(y~x, data=df[posList_ALL_Hypo,], pch=19, col="blue", cex=.2,type ="p")
  # text(y~x, data=df[posListLabels,], labels=dataPlotList[[i]][posListLabels,1], cex= 0.6,pos=2*(1+as.integer(df$x / df$y > 1)))
  
 # Create plot and capture as image and use to reduce reendering in the PDF
  require(grid)
  dev.new(noRStudioGD=T,height=7,width=7)
  par(mar=c(0,0,0,0),mai=c(0,0,0,0),pty="s")
  plot(y~x, data=df[order(df$dens),], col="dark grey",# col=cols, # Reordering rows so that densest points are plotted on top
       pch=19, cex=.2,type ="p",
       xlim=c(xmin,xmax),ylim=c(xmin,xmax),
       axes = F,frame.plot = F,asp=1,#bty ="n",
       xaxs = "i", yaxs = "i",xlab = "",ylab="")
  points(y~x, data=df[posList_ALL_Hyper,], pch=19, col="red", cex=.2,type ="p")
  points(y~x, data=df[posList_ALL_Hypo,], pch=19, col="blue", cex=.2,type ="p")
  capture <- grid.cap()
  dev.off()

  #Create PDF with plot + Frame
  PDFName <- paste("Results_RNABS/RNABS_2_Scatter_",ListNames[i],".pdf",sep="")
  pdf(PDFName,width=10,height=10,useDingbats=FALSE)
  par(pty="s")
  plot(1,type="n",
       xlim=c(xmin,xmax),ylim=c(xmin,xmax),
       axes = T,frame.plot = F,asp=1,
       ylab = paste(colnames(dataPlotList[[i]])[3],"(log2 RPM)") ,
       xlab = paste(colnames(dataPlotList[[i]])[2],"(log2 RPM)") )
  rasterImage(capture,xmin,xmin,xmax,xmax)
  text(y~x, data=df[posListLabels,], labels=dataPlotList[[i]][posListLabels,1], cex= 0.6,pos=2*(1+as.integer(df$x / df$y > 1)))
  legend(xmin,xmax, legend=c("Overlapping Hypermethylated DMRs", "Overlapping Hypomethylated DMRs"),
         col=c("red", "blue"), lty=1,lwd=2, cex=0.6,bty="n")
  #text(xmax,xmin,paste("R =",round(RSQ,3)),cex= 0.9,pos=2)
  
  dev.off()
}
