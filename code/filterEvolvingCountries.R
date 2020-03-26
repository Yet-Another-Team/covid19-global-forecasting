args = commandArgs(trailingOnly=TRUE)

require(stringr)

inDir <- args[1]
outDir <- args[2] # should end with *.csv

inFiles <- dir(path=inDir,pattern='*.csv$')

for(inFile in inFiles) {
  inFullPath <- file.path(inDir,inFile)
  outFullPath <- file.path(outDir,inFile)
  df <- read.csv(inFullPath)
  doCopy <- (length(df$Infected)>5) && (df$Infected[length(df$Infected)] > df$Infected[length(df$Infected)-5])
  if(doCopy) {
    print(paste0("Retain ",inFile,' as the epidemy still evolves at the location'))
    file.copy(inFullPath,outFullPath)
  }
}
print("Done")
