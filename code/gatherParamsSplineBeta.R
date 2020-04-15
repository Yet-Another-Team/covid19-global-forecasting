require(RJSONIO)
require(stringr)
args = commandArgs(trailingOnly=TRUE)

inputPattern <- args[1]
inputPattern <- substr(inputPattern,2,nchar(inputPattern)-1) # trimming out ''
outFile <- args[2]

asterixInPos = regexpr('\\*',inputPattern)[1]
print(paste0("asterisk pos is ",asterixInPos))
inDir <- substr(inputPattern,1,asterixInPos-2)
filemask <- substr(inputPattern,asterixInPos,nchar(inputPattern))

print(paste0("Looking for the files ",filemask," in ",inDir))
inFiles <- dir(path=inDir,pattern=paste0(filemask,'$'))

acc <- NULL

start_date <- as.Date(strptime('2020-01-01',format='%Y-%m-%d',tz="GMT"))

for(pars_df_path in inFiles) {
  print(paste0("Loading ",pars_df_path))
  inFilePath <- file.path(inDir,pars_df_path)
  json_data <- fromJSON(inFilePath)
  key <- substr(pars_df_path,1,nchar(pars_df_path)-5)
  atPos = regexpr('@',key)[1]
  if(atPos>=1) {
    province <- substr(key,1,atPos-1)
    country <- substr(key,atPos+1,nchar(key))
  } else {
    province <- ''
    country <- key
  }
  
  json_data <- within(json_data,rm(WeeklyBeta))
  json_data <- within(json_data,rm(WeeklyR0))
  json_data <- within(json_data,rm(BetaScales))
  json_data <- within(json_data,rm(RawGamma))
  json_data <- within(json_data,rm(bRaw))
  json_data <- within(json_data,rm(RawFirstDayInfectedCount))
  
  doSkip = F
  for(entry in json_data) {
    if(is.null(entry)) {
      print("null value in JSON Skippiing")
      doSkip <- T
      break
    }
  }
  if(doSkip)
    next
  
  
  #print(json_data)
  
  #pars_df <- data.frame(do.call("rbind",list(json_data)))
  pars_df <- as.data.frame(do.call("rbind",list(unlist(json_data))))
  
  origColsN <- ncol(pars_df)
  
  #print(pars_df)
  pars_df$Province <- province
  pars_df$Country <- country
  
  addedDateParams <- 0
  
  for(colName in names(pars_df)) {
    if(endsWith(colName,"DayNum")) {
      prefix <- substr(colName,0,nchar(colName) - nchar("DayNum"))
      outName <- paste0(prefix,"Date")
      dateVal <- as.integer(round(as.numeric(pars_df[[colName]])))
      dateVal <- as.character(start_date+dateVal-1) 
      pars_df[[outName]] <- dateVal
      addedDateParams <- addedDateParams + 1
    }
  }
  pars_df <- as.data.frame(pars_df)
  if(is.null(acc))
    acc <- pars_df
  else
    acc <- rbind(acc,pars_df)
  print(paste0(pars_df_path," processed"))
}
#acc<- as.data.frame(acc)
#print(acc)

write.csv(acc,file=outFile,row.names = F, quote = T)