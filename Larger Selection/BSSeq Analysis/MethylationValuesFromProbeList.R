filenameRAW <- file.choose() 
data.RAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")

data.RAW <- data.RAW[c((which (names(data.RAW) == "Distance")+1):length(data.RAW))]
se <- function(x) sqrt(var(x)/length(x))
data <-data.frame(mean=apply(data.RAW,2,mean),std=apply(data.RAW,2,sd),sem=apply(data.RAW,2,se))

