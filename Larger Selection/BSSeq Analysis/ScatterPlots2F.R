N <- data

#data <- data[1:100,]
i=4
j=3
for (i in 2:length(N)){
  for (j in 2:length(N)){
    pdf(paste(names(N)[i],"vs",names(N)[j],".pdf"),width=5,height=5,useDingbats=FALSE)
    par(pty="s")
  
    x1 <- N[[i]]
    x2 <- N[[j]]
    df <- data.frame(x1,x2)
    
    ## Use densCols() output to get density at each point
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    
    ## Map densities to colors
    cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                                "#FCFF00", "#FF9400", "#FF3100"))(256)
    df$col <- cols[df$dens]
    ## Plot it, reordering rows so that densest points are plotted on top
    plot(x2~x1, data=df[order(df$dens),], pch=19, col=col, cex=.2,type ="p",
         xlim=c(0,100),ylim=c(0,100),axes = T,frame.plot = F,asp=1,
         xlab = names(N)[i],ylab=names(N)[j],main="CGI : %CpG methylation")
    dev.off()
    
  }
}


