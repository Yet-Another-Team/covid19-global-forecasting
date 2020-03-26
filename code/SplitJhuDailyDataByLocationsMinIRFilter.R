args = commandArgs(trailingOnly=TRUE)

require(stringr)

jhuInputDir <- args[1]
populationCountFile <- args[2]
outputPattern <- args[3] # should end with *.csv
globalOutputFile <- args[4]
minAllowedConfirmedCount <- as.integer(args[5])
minAllowedRemovedCount <- as.integer(args[6])

asterixPos = regexpr('\\*',outputPattern)[1]
outDir <- substr(outputPattern,1,asterixPos-2)

population <- read.csv(populationCountFile)

dailyRepFiles <- dir(jhuInputDir,pattern='*.csv$')

sort(dailyRepFiles)

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
  names(dailyDfAgg)[1:2] <- c('Province_State','Country_Region')
  
  dateName <- substr(dailyRepFile,1,nchar(dailyRepFile)-4)
  start_date <- strptime('01-01-2020',format='%m-%d-%Y',tz="GMT")
  cur_date <- strptime(dateName,format='%m-%d-%Y',tz="GMT")
  cur_day <- difftime(cur_date,start_date,units='days')+1
  
  dailyDfAgg$date <- cur_date
  dailyDfAgg$dayNum <- cur_day
  
  N <- nrow(dailyDfAgg)
  for(i in (1:N)) {
    row <- dailyDfAgg[i,]
    province <- row$Province_State
    country <- row$Country_Region
    key <- paste0(province,'@',country)
    curValLatestUpdate <- row$Last.Update
    #print(curValLatestUpdate)
    #print(row)
    keyDf <- locationList[[key]]
    #latestUpdateAtLoc <- locationLatest[[key]]
    
    # whether to register this day or not
    #validEntry <- (is.null(latestUpdateAtLoc)) || (difftime(curValLatestUpdate,latestUpdateAtLoc,units='days')>0)
    validEntry <- TRUE
    
    if(validEntry) {
      if(is.null(keyDf)) {
        keyDf <- row
      }
      else
        keyDf <- rbind(keyDf,row)
      locationList[[key]] <- keyDf
      #latestUpdateAtLoc[[key]] <- curValLatestUpdate
    } else {
      print(paste0("skiping day ",cur_date," for ",key," as there were no updates"))
    }
  }
}

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

for(i in (1:nrow(population))) {
  province <- as.character(population[i,1])
  country <- as.character(population[i,2])
  popCount <- population[i,3]
  #print(paste0("Extracting province '",province,"' state '",country,"' with population ",popCount))
  
  siteDf <- get_site_df(province,country,popCount)
  
  #print(head(siteDf))
  if(is.null(siteDf)) {
    print(paste0("Can't find data for province '",province,"' country '",country,"'"))
    next;
  }
  
  if(is.null(globalTs))
    globalTs <- siteDf[,3:10]
  else
    globalTs <- rbind(globalTs,siteDf[,3:10])
  
  maxConfirmed <- max(siteDf$Confirmed)
  #if(is.na(maxConfirmed)) {
  #  print(paste0(province,' - ',country))
  #  print(df$Confirmed)
  #}
  maxRemoved <- max(siteDf$Removed)
  if(maxConfirmed < minAllowedConfirmedCount) {
    print(paste0("Skipping province '",province,"' country '",country,"' due to too low confirmed cases: ",maxConfirmed))
    next; 
  }
  if(maxRemoved < minAllowedRemovedCount) {
    print(paste0("Skipping province '",province,"' country '",country,"' due to too low removed cases: ",maxRemoved))
    next; 
  }
  if(nchar(province)>0)
    province <- paste0(province,'@')
  key <- paste0(province, country)
  key <- str_replace_all(key," ","_")
  key <- str_replace_all(key,"\\*","#")
  outFile <- file.path(outDir,paste0(key,'.csv'))
  #print(paste0('Writing ',outFile))
  write.csv(siteDf,outFile, row.names = F)
}
#print(head())

globalTsAgg <- aggregate(cbind(Confirmed,Active,Infected,Deaths,Recovered,Removed) ~ dayNum,data = globalTs, sum)

#globalTsAgg <- globalTsAgg[globalTsAgg$dayNum >=71,] #removing due to QC

#writing global
globalTsAgg$Susceptible <- globalPopulation - globalTsAgg$Infected - globalTsAgg$Removed
write.csv(globalTsAgg,globalOutputFile, row.names = F)

print("Done")
