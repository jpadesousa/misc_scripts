#biocLite("BSgenome.Mmusculus.UCSC.mm10")


library(BSgenome.Mmusculus.UCSC.mm10)

chromNames <- names(Mmusculus)[!grepl("_", names(Mmusculus))]

CGs <- lapply(chromNames, function(x) start(matchPattern("CG", Mmusculus[[x]])))

names(CGs) <- chromNames

str(CGs)

start(matchPattern("CG", Mmusculus[[1]])) -> b
