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

provinces <- c()
countries <- c()
deaths <- c()
recovered <- c()

print(paste0('Will process ',length(inFiles),' files'))

for(pars_df_path in inFiles) {
  inFilePath <- file.path(inDir,pars_df_path)
  df1 <- read.csv(inFilePath)
  N <- nrow(df1)
  
  key <- substr(pars_df_path,1,nchar(pars_df_path)-4)
  atPos = regexpr('@',key)[1]
  if(atPos>=1) {
    province <- substr(key,1,atPos-1)
    country <- substr(key,atPos+1,nchar(key))
  } else {
    province <- ''
    country <- key
  }
  
  print(paste0('file ',pars_df_path,' province ',province,' country ',country))
  
  provinces <- c(provinces, province)
  countries <- c(countries, country)
  
  deaths <- c(deaths, df1$Deaths[N])
  recovered <- c(recovered, df1$Recovered[N])
}



resDf <- data.frame(Province=provinces, Country=countries, Deaths = deaths, Recovered = recovered)

print(paste0(sum(resDf$Deaths>0),' locations are to be discarded due to zero deaths:'))
discardedDf <- resDf[resDf$Deaths>0,1:2]
print(discardedDf)

resDf <- resDf[,]
resDf$DeathPortion <- resDf$Deaths / (resDf$Deaths +resDf$Recovered)

print(resDf)

write.csv(resDf, file=outFile, row.names = F)
