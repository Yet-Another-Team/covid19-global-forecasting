args = commandArgs(trailingOnly=TRUE)

require(stringr)
require(pracma)

jhuInputDir <- args[1]
populationCountFile <- args[2]
inputKeyFilterFile <- args[3]
outputPattern <- args[4] # should end with *.csv
globalOutputFile <- args[5]
minAllowedConfirmedCount <- as.integer(args[6])
minAllowedRemovedCount <- as.integer(args[7])

asterixPos = regexpr('\\*',outputPattern)[1]
outDir <- substr(outputPattern,1,asterixPos-2)

population <- read.csv(populationCountFile)
population$Province.State <- str_replace_all(population$Province.State," ","_")
population$Country.Region <- str_replace_all(population$Country.Region," ","_")
population$Province.State <- str_replace_all(population$Province.State,"\\*","#")
population$Country.Region <- str_replace_all(population$Country.Region,"\\*","#")


dailyRepFiles <- dir(jhuInputDir,pattern='*.csv$')

sort(dailyRepFiles)

inputKeyFilterDf <- read.csv(inputKeyFilterFile)
neededKeys <- inputKeyFilterDf$Key
  

locationList <- list()
#locationLatest <- list()

for(dailyRepFile in dailyRepFiles) {
  dailyDf <- read.csv(file.path(jhuInputDir,dailyRepFile), fileEncoding = 'UTF-8-BOM')
  if(sum(is.na(dailyDf$Confirmed))>0)
    dailyDf[is.na(dailyDf$Confirmed),]$Confirmed <- 0
  if(sum(is.na(dailyDf$Deaths))>0)
    dailyDf[is.na(dailyDf$Deaths),]$Deaths <- 0
  if(sum(is.na(dailyDf$Recovered))>0)
    dailyDf[is.na(dailyDf$Recovered),]$Recovered <- 0
  if('Province.State' %in% names(dailyDf)) {
    dailyDfAgg <- aggregate(cbind(Confirmed,Deaths,Recovered) ~ Province.State + Country.Region,data = dailyDf, sum)
  } else {
    dailyDfAgg <- aggregate(cbind(Confirmed,Deaths,Recovered,Active) ~ Province_State + Country_Region,data = dailyDf, sum)
  }
  dailyDfAgg$Removed <- dailyDfAgg$Deaths + dailyDfAgg$Recovered
  dailyDfAgg$Active <- dailyDfAgg$Confirmed - dailyDfAgg$Removed
  dailyDfAgg$Infected <- dailyDfAgg$Active
  names(dailyDfAgg)[1:2] <- c('Province.State','Country.Region')
  
  dateName <- substr(dailyRepFile,1,nchar(dailyRepFile)-4)
  start_date <- strptime('01-01-2020',format='%m-%d-%Y',tz="GMT")
  cur_date <- strptime(dateName,format='%m-%d-%Y',tz="GMT")
  cur_day <- difftime(cur_date,start_date,units='days')+1
  
  dailyDfAgg$date <- cur_date
  dailyDfAgg$dayNum <- cur_day
  
  N <- nrow(dailyDfAgg)
  for(i in (1:N)) {
    row <- dailyDfAgg[i,]
    province <- as.character(row$Province.State)
    country <- as.character(row$Country.Region)
    if(nchar(province)>0)
      province <- paste0(province,'@')
    key <- paste0(province, country)
    key <- str_replace_all(key," ","_")
    key <- str_replace_all(key,"\\*","#")
    
    curValLatestUpdate <- row$Last.Update
    #print(curValLatestUpdate)
    #print(row)
    keyDf <- locationList[[key]]
    #latestUpdateAtLoc <- locationLatest[[key]]
    
    
    if(is.null(keyDf)) {
      keyDf <- row
    }
    else
      keyDf <- rbind(keyDf,row)
    locationList[[key]] <- keyDf
    #latestUpdateAtLoc[[key]] <- curValLatestUpdate
  
  }
}

print("Loaded all daily reports files")

print(paste0("There are ",length(locationList),' locations in data location list'))
print(paste0("There are ",nrow(population),' locations in location population df'))

# N is the population of the country
get_site_df <- function(province,country, N) {
  key <- paste0(province,'@',country)
  df <- locationList[[key]]
  if(!is.null(df))
    df$Susceptible <- N - df$Infected - df$Removed
  return(df)
}

globalTs <- NULL
globalPopulation <- 7700000000

#print(population[,1:2])
for(key in neededKeys) {
  atPos = regexpr('\\@',key)[1]
  if(atPos>0) {
    province <- as.character(substr(key, 1, atPos-1))
    country <- as.character(substr(key,atPos+1,nchar(key)))
  } else {
    province <- ''
    country <- key
  }
  popRow <- population[(population$Province.State == province) & (population$Country.Region == country),]
  if(nrow(popRow)>0) {
    popCount <- population[1,3]
    siteDf <- get_site_df(province,country,popCount)
    
    #print(head(siteDf))
    if(is.null(siteDf)) {
      print(paste0("!!! ---> Can't find data for province '",province,"' country '",country,"'"))
      next;
    }
    
    if(is.null(globalTs))
      globalTs <- siteDf[,3:10]
    else
      globalTs <- rbind(globalTs,siteDf[,3:10])
    
    maxConfirmed <- max(siteDf$Confirmed)
    maxRemoved <- max(siteDf$Removed)
    if(maxConfirmed < minAllowedConfirmedCount) {
      print(paste0("Skipping province '",province,"' country '",country,"' due to too low confirmed cases: ",maxConfirmed))
      next; 
    }
    if(maxRemoved < minAllowedRemovedCount) {
      print(paste0("Skipping province '",province,"' country '",country,"' due to too low removed cases: ",maxRemoved))
      next; 
    }
    infMax <- max(siteDf$Infected)
    #print(df$Infected/infMax*10.0)
    infCoarse <- round(siteDf$Infected/infMax*30.0)
    infectedMA <- movavg(infCoarse, 10, type="t")
    N <- length(infectedMA)
    if(N<3) {
      print(paste0("Skipping province '",province,"' country '",country,"' due to too few data points: ",N))
      next;
    }
    infectedMaDeriv <- infectedMA[2:N] - infectedMA[1:(N-1)]
    extremums <- sign(infectedMaDeriv[2:(N-1)]*infectedMaDeriv[1:(N-2)]) == -1
    extremumCount <- sum(extremums) # derivative changes its sign
    
    if(extremumCount>1) {
      print(paste0("-----> Skipping province '",province,"' country '",country,"' due to too many infected extremum count: ",extremumCount," (e.g. multiple peaks)"))
      # print(extremums)
      next;
    }
    else {
      #print(paste0("province '",province,"' country '",country,"' extremum count: ",extremumCount))
      #print(infectedMaDeriv)
    }
    
    outFile <- file.path(outDir,paste0(key,'.csv'))
    #print(paste0('Writing ',outFile))
    write.csv(siteDf,outFile, row.names = F)
    print(paste0("Done with province '",province,"' state '",country,"' with population ",popCount))
  } else {
    print(paste0("!!! ---> Can't find population data for '",province,"' state of '",country,"'"))
  }
}
#print(head())

globalTsAgg <- aggregate(cbind(Confirmed,Active,Infected,Deaths,Recovered,Removed) ~ dayNum,data = globalTs, sum)

#globalTsAgg <- globalTsAgg[globalTsAgg$dayNum >=71,] #removing due to QC

#writing global
globalTsAgg$Susceptible <- globalPopulation - globalTsAgg$Infected - globalTsAgg$Removed
write.csv(globalTsAgg,globalOutputFile, row.names = F)

print("Done")
