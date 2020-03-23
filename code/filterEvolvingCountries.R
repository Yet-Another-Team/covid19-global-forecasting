args = commandArgs(trailingOnly=TRUE)

require(stringr)

inDir <- args[1]
outDir <- args[2] # should end with *.csv

inFiles <- dir(path=inDir,pattern='*.csv$')

for(inFile in inFiles) {
  inFullPath <- file.path(inDir,inFile)
  outFullPath <- file.path(outDir,inFile)
  df <- read.csv(inFullPath)
  doCopy <- df$infected[length(df$infected)] > median(df$infected)
  if(doCopy) {
    print(paste0("Retain ",inFile,' as the epidemy still evolves at the location'))
    file.copy(inFullPath,outFullPath)
  }
}
print("Done")
