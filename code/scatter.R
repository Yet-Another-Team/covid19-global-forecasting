args = commandArgs(trailingOnly=TRUE)

inDir <- args[1]
outputPattern <- args[2] # should end with *.csv
outputPattern <- substr(outputPattern,2,nchar(outputPattern)-1) # trimming out ''

print(paste0("outputPattern is ",outputPattern))
print(paste0("inputDir is ",inDir))

asterixOutPos = regexpr('\\*',outputPattern)[1]

print(paste0("asterisk pos is ",asterixOutPos))
outDir <- substr(outputPattern,1,asterixOutPos-2)
filemask <- substr(outputPattern,asterixOutPos,nchar(outputPattern))

print(paste0("Generating files ",filemask," in ",outDir))
inFiles <- dir(path=inDir,pattern=paste0(filemask,'$'))
print(paste0(length(inFiles)," to copy"))
for(inFile in inFiles) {
  inFullPath <- file.path(inDir,inFile)
  outFullPath <- file.path(outDir,inFile)
  print(paste0("Copiyng ",inFile))
  file.copy(inFullPath,outFullPath)
}
print("Done")
