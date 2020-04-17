args = commandArgs(trailingOnly=TRUE)

require(RJSONIO)
require(stringr)

outputPattern <- args[1] # should end with *.csv
outputPattern <- substr(outputPattern,2,nchar(outputPattern)-1) # trimming out ''
#print(paste0("outputPattern is ",outputPattern))

asterixOutPos = regexpr('\\*',outputPattern)[1]

#print(paste0("asterisk pos is ",asterixOutPos))
outDir <- substr(outputPattern,1,asterixOutPos-2)
filemask <- substr(outputPattern,asterixOutPos,nchar(outputPattern))
print(paste0("out dir: ",outDir," ; filemask: ",filemask))

locationsKeys <- c(
  'Austria',
  'Beijing@China',
  'Heilongjiang@China',
  'Italy',
  'Korea,_South'
)
sigmaValues <- c(12,4,2,1,0.5, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8,1/10,1/14,1/20)

print(paste0("Locations to test: ",length(locationsKeys)))
print(paste0("Sigma values to test: ",length(sigmaValues)))
print(paste0("Total experiments count: ",length(sigmaValues) * length(locationsKeys)))
for(locKey in locationsKeys)
  for(sigma in sigmaValues) {
    outFile <- file.path(outDir,paste0(locKey,'-',str_replace_all(as.character(round(sigma,3)),"[.]","_"),'.json'))
    outList <- list(
      locationKey = locKey,
      sigma = sigma
    )
    exportJson <- toJSON(outList)
    write(exportJson, outFile)
  }
print("Done")