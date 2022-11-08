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
i=1


#install.packages("gProfileR")
library("gProfileR")

for (i in 1:4) {
  
  print (ListNames[i])
  # DMRGeneListHyper
  OutputName <- paste ("Results_RNABS/HyperDMR_Exp_Gene_",ListNames[i],".txt",sep="")
  Selection <- (match(ToGOListALL[[i]],DMRGeneListHyper) > 0)
  
  Output <- na.omit(dataset[(match(rownames(dataset),na.omit(ToGOListALL[[i]][Selection])) > 0),])
  Output <- cbind('Gene' = rownames(Output),
                  'log2 Fold PCE vs Ctrl (CE)' = apply(Output[,c(1:6)],1,mean) - apply(Output[,c(7:12)],1,mean),
                  'log2 Fold PCE vs Ctrl (RT)' = apply(Output[,c(13:18)],1,mean) - apply(Output[,c(19:24)],1,mean),  
                  'BAT: PCE-CE (Avg)' =  apply(Output[,c(1:6)],1,mean),
                  'BAT: Ctrl-CE (Avg)' = apply(Output[,c(7:12)],1,mean),
                  'BAT: PCE-RT (Avg)' =  apply(Output[,c(13:18)],1,mean),
                  'BAT: Ctrl-RT (Avg)' =  apply(Output[,c(19:24)],1,mean))
  
  
  write.table(Output , OutputName, sep="\t",row.names = FALSE)


  
  
  # DMRGeneListHyo
  Selection <- (match(ToGOListALL[[i]],DMRGeneListHypo) > 0 )
  OutputName <- paste ("Results_RNABS/HypoDMR_Exp_Gene_",ListNames[i],".txt",sep="")

  Output <- na.omit(dataset[(match(rownames(dataset),na.omit(ToGOListALL[[i]][Selection])) > 0),])
  Output <- cbind('Gene' = rownames(Output),
                  'log2 Fold PCE vs Ctrl (CE)' = apply(Output[,c(1:6)],1,mean) - apply(Output[,c(7:12)],1,mean),
                  'log2 Fold PCE vs Ctrl (RT)' = apply(Output[,c(13:18)],1,mean) - apply(Output[,c(19:24)],1,mean),  
                  'BAT: PCE-CE (Avg)' =  apply(Output[,c(1:6)],1,mean),
                  'BAT: Ctrl-CE (Avg)' = apply(Output[,c(7:12)],1,mean),
                  'BAT: PCE-RT (Avg)' =  apply(Output[,c(13:18)],1,mean),
                  'BAT: Ctrl-RT (Avg)' =  apply(Output[,c(19:24)],1,mean))
  
  
  write.table(Output , OutputName, sep="\t",row.names = FALSE)
  
  
  
  # DMRGeneListHyper
  # GO_Hyper <- gprofiler(ToGOListALL[[i]][SelectionHyper],organism = "mmusculus",png_fn = NULL, 
  #                     src_filter = c("GO:BP", "GO:MF", "GO:CC","KEGG"))
  # OutputName <- paste ("Results_RNABS/HyperDMR_Exp_GO_",ListNames[i],".txt",sep="")
  # write.table(GO_Hyper , OutputName, sep="\t")
  # 
  # # DMRGeneListHypo
  # GO_Hypo <- gprofiler(ToGOListALL[[i]][SelectionHypo],organism = "mmusculus",png_fn = NULL, 
  #                     src_filter = c("GO:BP", "GO:MF", "GO:CC","KEGG"))
  # OutputName <- paste ("Results_RNABS/HypoDMR_Exp_GO_",ListNames[i],".txt",sep="")
  # write.table(GO_Hypo, OutputName, sep="\t")
  
  
}
