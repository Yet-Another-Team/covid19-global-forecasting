# initial sigmoid 
# y = sigmoid((days-inflectionDayNum)/k)*EstimatedSusceptiblePopulation
# y = sigmoid(x) * EstimatedSusceptiblePopulation;
# x = (days-inflectionDayNum)/k

 
# back transform is
# x_back = x*k + inflectionDayNum

require(stringr)

args = commandArgs(trailingOnly=TRUE)

sigmoidParamsFile <- args[1]
outputFile <- args[2]
argc <- length(args)
tsSourceDirs<-args[3:argc]

print("TS source dirs:")
print(tsSourceDirs)

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

outDf <- NULL

sigmoidParamsDf <- read.csv(sigmoidParamsFile)
N <- nrow(sigmoidParamsDf)

for(i in (1:N)) {
  dfRow <- sigmoidParamsDf[i,]
  inflDayNum <- dfRow$ConfirmedInflectionDayNum
  k <- dfRow$ConfirmedK
  susceptPop <- dfRow$EstimatedSusceptiblePopulation
  
  province <- dfRow$Province
  country <- dfRow$Country
  
  key <- getKey(province,country)
  tsFileName <- paste0(key,'.csv')
  
  tsDays <- NULL
  tsDeaths <- NULL
  for(tsSourceDir in tsSourceDirs) {
    tsFileName2 <- file.path(tsSourceDir,tsFileName)
    #print(paste0('Tring to find ',tsFileName2))
    if(file.exists(tsFileName2)) {
      tsDf <- read.csv(tsFileName2)
      tsDays <- tsDf$dayNum
      tsDeaths <- tsDf$Deaths
      break;
    }
  }
  if(is.null(tsDays)) {
    print(paste0("WARNING: could not find timeseries for ",key))
    stop()
  }
  
  # back transform is
  # y_back = y/EstimatedSusceptiblePopulation
  # x_back = x*k + inflectionDayNum
  sigmoidSpaceY <- tsDeaths/susceptPop
  sigmoidSpaceX <- (tsDays - inflDayNum)/k
  
  curLocDf <- data.frame(x = sigmoidSpaceX, y = sigmoidSpaceY, province=province, country=country )
  print(paste0("Got ",nrow(curLocDf)," points for ",key))
  if(is.null(outDf))
    outDf <- curLocDf
  else
    outDf <- rbind(outDf,curLocDf)
}
write.csv(outDf,file=outputFile,row.names = F)
print("Done")
