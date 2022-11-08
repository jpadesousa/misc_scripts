
data.plot <- na.omit(data[,c(1,4,7,8)])

,"WT Genebody",BSLab = "Genebody methylation", PDF=T,Reg=F)


data.plot <- data.plot[data.plot[2]>min(data.plot[2]),]
data.plot <- data.plot[data.plot[3]!=0,]
names(data.plot) <- c("Gene","RNA","BS","CGI")

x1 <- data.plot$BS
x2 <- data.plot$RNA
df <- data.frame(x1,x2)

## Use densCols() output to get density at each point
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
x <- 1+ as.numeric(data.plot$CGI)
df$col <- x

#Plot
plot(x2~x1, data=df[order(data.plot$CGI),], col=col,pch=19, cex=.2,type ="p",
     axes = T,frame.plot = F,
     ylab="rel. mRNA expression")
