library(epiGTBN)
#========================================================
# load parameters for GTBN
source("parameters.R")
#========================================================
# read data
mydata<- read.table("0.4HER0.4MAF3LOC_EDM-1_004.txt", header=TRUE)
#========================================================
# change data type from integer to numeric (e.g change data's class from integer to double)
tmp <- data.frame(matrix(as.numeric(unlist(mydata)), ncol = length(mydata[1,])))
rownames(tmp) <- rownames(mydata)
colnames(tmp) <- colnames(mydata)
mydata <- tmp
#========================================================
# if necessary, data can be discretized using discretize()
# mydata <- discretize(mydata, method = 'interval', breaks=3)
#========================================================
timestart <- Sys.time()

# sink(file = "GTBN-3nodes-log.txt")
res <- gtbn3(mydata, max.iter = 60, debug = FALSE)
# sink()

timeend <- Sys.time()
time_gtbn <- c(timeend-timestart)
#========================================================
prediction <- data.frame(res$arcs)