# #########################################################  
# GO Analysis using gProfileR
# #########################################################  


# Get Gene Lists

setwd("/Users/meyennf/Cambridge/Projects/Collaboration_Wolfrum/2017_PCE/2017October/")

GeneListFolderMed <- "Results_RNA/7_RNA_Lists/2-MedStringency/"
GeneListFileNamesMed <- list.files(GeneListFolderMed, pattern=".txt")
ListNames <- unlist(strsplit(GeneListFileNamesMed,"_Med.txt"))

GeneList <- list()
ToGOListALL <- list()
ToGOListUP <- list()
ToGOListDOWN <- list()

for (i in 1:4) {
  print(ListNames[i])
  GeneList[[ListNames[i]]] <- read.csv(paste(GeneListFolderMed,GeneListFileNamesMed[i],sep=""),
                                          sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")
  ToGOListALL[[ListNames[i]]] <- as.character(GeneList[[i]][,1]) 
  }

# 1 Ctrl: 3CE vs RT: CTRL_3CEvsRT
ToGOListUP$CTRL_3CEvsRT   <-  as.character(GeneList$CTRL_3CEvsRT[(GeneList$CTRL_3CEvsRT$CTRL_3CE - GeneList$CTRL_3CEvsRT$CTRL_RT > 0),1])
ToGOListDOWN$CTRL_3CEvsRT <-  as.character(GeneList$CTRL_3CEvsRT[(GeneList$CTRL_3CEvsRT$CTRL_3CE - GeneList$CTRL_3CEvsRT$CTRL_RT < 0),1])

# 2 PCE: 3CE vs RT: PCE_3CEvsRT
ToGOListUP$PCE_3CEvsRT   <-  as.character(GeneList$PCE_3CEvsRT[(GeneList$PCE_3CEvsRT$PCE_3CE - GeneList$PCE_3CEvsRT$PCE_RT > 0),1])
ToGOListDOWN$PCE_3CEvsRT <-  as.character(GeneList$PCE_3CEvsRT[(GeneList$PCE_3CEvsRT$PCE_3CE - GeneList$PCE_3CEvsRT$PCE_RT < 0),1])

# 3 PCE vs CTRL (at 3CE): PCEvsCTRL_3CE
ToGOListUP$PCEvsCTRL_3CE   <-  as.character(GeneList$PCEvsCTRL_3CE[(GeneList$PCEvsCTRL_3CE$PCE_3CE - GeneList$PCEvsCTRL_3CE$CTRL_3CE > 0),1])
ToGOListDOWN$PCEvsCTRL_3CE <-  as.character(GeneList$PCEvsCTRL_3CE[(GeneList$PCEvsCTRL_3CE$PCE_3CE - GeneList$PCEvsCTRL_3CE$CTRL_3CE < 0),1])

# 4 PCE vs CTRL (at RT): PCEvsCTRL_RT
ToGOListUP$PCEvsCTRL_RT   <-  as.character(GeneList$PCEvsCTRL_RT[(GeneList$PCEvsCTRL_RT$PCE_RT - GeneList$PCEvsCTRL_RT$CTRL_RT > 0),1])
ToGOListDOWN$PCEvsCTRL_RT <-  as.character(GeneList$PCEvsCTRL_RT[(GeneList$PCEvsCTRL_RT$PCE_RT - GeneList$PCEvsCTRL_RT$CTRL_RT < 0),1])



# GO Analysis

#install.packages("gProfileR")
library("gProfileR")

for (i in 1:4) {
  
  print (ListNames[i])
  

  # Combined Up and Down Expression
  GO_All <- gprofiler(ToGOListALL[[i]],organism = "mmusculus",png_fn = NULL, 
                  src_filter = c("GO:BP", "GO:MF", "GO:CC","KEGG"))
  OutputName <- paste ("Results_RNA/6_RNA_GO/2-MedStringency/GO_",ListNames[i],"_both.txt",sep="")
  write.table(GO_All, OutputName, sep="\t")
  
  # Up Expression
  GO_Up <- gprofiler(ToGOListUP[[i]],organism = "mmusculus",png_fn = NULL, 
                     src_filter = c("GO:BP", "GO:MF", "GO:CC","KEGG"))
  OutputName <- paste ("Results_RNA/6_RNA_GO/2-MedStringency/GO_",ListNames[i],"_UP.txt",sep="")
  write.table(GO_Up, OutputName, sep="\t")
  
  # Down Expression  
  GO_Down <- gprofiler(ToGOListDOWN[[i]],organism = "mmusculus",png_fn = NULL, 
                       src_filter = c("GO:BP", "GO:MF", "GO:CC","KEGG"))
  OutputName <- paste ("Results_RNA/6_RNA_GO/2-MedStringency/GO_",ListNames[i],"_DOWN.txt",sep="")
  write.table(GO_Down, OutputName, sep="\t")
  
  
}







