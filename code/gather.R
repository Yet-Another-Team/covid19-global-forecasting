args = commandArgs(trailingOnly=TRUE)

inputPattern <- args[1]
inputPattern <- substr(inputPattern,2,nchar(inputPattern)-1) # trimming out ''
outDir <- args[2] # should end with *.csv

print(paste0("input pattern is ",inputPattern))
print(paste0("output dir is ",outDir))

asterixInPos = regexpr('\\*',inputPattern)[1]
print(paste0("asterisk pos is ",asterixInPos))
inDir <- substr(inputPattern,1,asterixInPos-2)
filemask <- substr(inputPattern,asterixInPos,nchar(inputPattern))

print(paste0("Looking for the files ",filemask," in ",inDir))
inFiles <- dir(path=inDir,pattern=paste0(filemask,'$'))
print(paste0(length(inFiles)," to copy"))
for(inFile in inFiles) {
  inFullPath <- file.path(inDir,inFile)
  outFullPath <- file.path(outDir,inFile)
  print(paste0("Copiyng ",inFile))
  file.copy(inFullPath,outFullPath,overwrite=T)
}
print("Done")
