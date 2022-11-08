require(grid)

require(Rtsne)
require(tsne)

require(calibrate)

setwd("~/Cambridge/Projects/Collaborations/Christian Wolfrum/PCE/Methylation/")
filenameRAW <- "AllValues_3&3.txt"

dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".")[,c(14:25)]

Rtsne(t(dataRAW), perplexity=3, max_iter=1000, verbose=T,PCA=T,check_duplicates = T,theta=0) -> tsne.results
plot(tsne.results$Y)
textxy(tsne.results$Y[1:6,1],tsne.results$Y[1:6,2],names(dataRAW)[1:6],col="red")
textxy(tsne.results$Y[7:12,1],tsne.results$Y[7:12,2],names(dataRAW)[7:12],col="blue")


cols <- rep(c("blue","red"),each=6)

# this is the epoch callback function used by tsne. 
# x is an NxK table where N is the number of data rows passed to tsne, and K is the dimension of the map. 
# Here, K is 2, since we use tsne to map the rows to a 2D representation (map).

#ecb = function(x, y){ plot(x, t='n'); text(x, labels=names(dataRAW), col=cols); }

#tsne_5MC = tsne(t(dataRAW[1:1000,]), k=2, max_iter=1000, epoch_callback = ecb, epoch=200, perplexity=50)




?tsne()


proc.time()
