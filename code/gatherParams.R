require(RJSONIO)
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
  
  pars_df <- data.frame(do.call("rbind",list(json_data)))
  origColsN <- ncol(pars_df)
  pars_df$FirstDate <- as.character(start_date+pars_df$FirstDayNum-1)
  pars_df$PeakDate <- as.character(start_date+pars_df$PeakDayNum-1)
  pars_df$Province <- province
  pars_df$Country <- country
  pars_df <- pars_df[,c(origColsN+3,origColsN+4,1:4,origColsN+1,5:6,origColsN+2,7:origColsN)]
  if(is.null(acc))
    acc <- pars_df
  else
    acc <- rbind(acc,pars_df)
}

write.csv(acc,file=outFile,row.names = F, quote = T)