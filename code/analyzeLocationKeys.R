require(stringr)

args = commandArgs(trailingOnly=TRUE)

jhuDrDir <- args[1]
jhuTsFile <- args[2]
kaggleTestFile <- args[3]
populationFile <-args[4]
neededFromDrFile <- args[5]
outDir <- args[6]

drFiles <- dir(jhuDrDir)
drFiles <- drFiles[endsWith(drFiles,".csv")]
drDf <- data.frame(files = drFiles)
drDf$dateStr <- substr(drDf$files,1,10)
drDf$date <- strptime(drDf$dateStr,format='%m-%d-%Y',tz="GMT")
maxDate <- max(drDf$date)
maxFile <- as.character(drDf$files[drDf$date == maxDate])

getKey <- function(province,country) {
  province <- as.character(province)
  country <- as.character(country)
  nonEmptyProvince <- nchar(province)>0
  province[nonEmptyProvince] <- paste0(province[nonEmptyProvince],'@')
  key <- paste0(province, country)
  key <- str_replace_all(key," ","_")
  key <- str_replace_all(key,"\\*","#")
  key
}


drDf <- read.csv(file.path(jhuDrDir,maxFile))
drKeys <- unique(getKey(drDf$Province_State, drDf$Country_Region))

print(paste0("Daily reports locations: ",length(drKeys)))

tsDf <- read.csv(jhuTsFile)
tsKeys <- unique(getKey(tsDf$Province.State,tsDf$Country.Region))
print(paste0("Time series reports locations: ",length(tsKeys)))

kaggleDf <- read.csv(kaggleTestFile)
kaggleKeys <- unique(getKey(kaggleDf$Province_State, kaggleDf$Country_Region))
print(paste0("Kaggle challange locations: ",length(kaggleKeys)))

popDf <- read.csv(populationFile)
popKeys <- unique(getKey(popDf$Province.State, popDf$Country.Region))
print(paste0("Population locations: ",length(popKeys)))

neededFromDr <- read.csv(neededFromDrFile)
neededFromDrKeys <- neededFromDr[,1]
print(paste0("Needed from dr locations:", length(neededFromDrKeys)))


missingPopKeys <- setdiff(kaggleKeys,popKeys)
print("")
print("Following keys are missing in population file:")
print("")
print(missingPopKeys)
write.csv(missingPopKeys,file=file.path(outDir,"missingPopKeys.csv"), row.names = F)

notCoveredKaggleKeys <- setdiff(kaggleKeys,tsKeys)
notCoveredKaggleKeys <- setdiff(notCoveredKaggleKeys, neededFromDrKeys)
print("")
print("Following Kaggle keys are not covered:")
print("")
print(notCoveredKaggleKeys)
write.csv(missingPopKeys,file=file.path(outDir,"notCoveredKaggleKeys.csv"), row.names = F)
print("Done")
