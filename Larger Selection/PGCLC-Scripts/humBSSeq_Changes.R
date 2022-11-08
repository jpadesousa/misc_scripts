
setwd("~/Documents/Cambridge/Projects/PGCLC/hBSSeq/")

filenameRAW <- "50CpG_5mC_Values5_3.txt"
dataRAW <- read.csv(filenameRAW,sep="\t",header=TRUE,skip=0,fill=TRUE,dec = ".",
                    as.is=T)
                    nrows = 10000)

data <- dataRAW[,c(2,3,4,13:length(dataRAW))]

data <- dataRAW[,c(2,3,4,14,16:21)]

# names(data)

# Remove Rows with NAs
data <- na.omit(data)
# Remove Zero Cols for  "H9.reset..merge."
data <- data [data$H9.reset..merge. != 0,]

orig <- data[,c(4:length(data))]
divided <- data[,c(4:length(data))] / data$H9.reset..merge.
substr <- data[,c(4:length(data))] - data$H9.reset..merge.

data2Use <- divided

# Remove Rows with NAs
data2Use <- na.omit(data2Use)
# Remove Zero Rows
data2Use <- data2Use[apply(data2Use[,2:length(data2Use)],1,sum) !=0 ,]



#Random selection of probes
sampling <- 10000
if (nrow(data2Use)< 10000) { sampling <- nrow(data2Use)}
data2Use <- data2Use[sample(nrow(data2Use), sampling), ]

data_matrix <- data.matrix(data2Use[,2:length(data2Use)])
