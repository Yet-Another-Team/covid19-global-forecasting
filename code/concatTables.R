args = commandArgs(trailingOnly=TRUE)

outFile <- args[1]
acc <- NULL
N <- length(args)
for(i in (2:N)) {
  inFilePath <- args[i]
  print(paste0("Reading table ",inFilePath))
  curDf <- read.csv(inFilePath)
  if(is.null(acc))
    acc <- curDf
  else
    acc <- rbind(acc,curDf)
}
print(paste0("Writing rbinded results into ",outFile))
write.csv(acc,file=outFile, row.names = F)