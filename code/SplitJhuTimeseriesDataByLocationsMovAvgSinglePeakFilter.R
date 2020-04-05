args = commandArgs(trailingOnly=TRUE)

require(stringr)
require(pracma)

jhuInputDir <- args[1]
populationCountFile <- args[2]
outputPattern <- args[3] # should end with *.csv
globalOutputFile <- args[4]
minAllowedConfirmedCount <- as.integer(args[5])
minAllowedRemovedCount <- as.integer(args[6])

asterixPos = regexpr('\\*',outputPattern)[1]
outDir <- substr(outputPattern,1,asterixPos-2)

jhu_confirmed <- read.csv(file.path(jhuInputDir,'time_series_covid19_confirmed_global.csv'),fileEncoding = "UTF-8-BOM")
jhu_deaths <- read.csv(file.path(jhuInputDir,'time_series_covid19_deaths_global.csv'),fileEncoding = "UTF-8-BOM")
jhu_recovered <- read.csv(file.path(jhuInputDir,'time_series_covid19_recovered_global.csv'),fileEncoding = "UTF-8-BOM")

#print(jhu_recovered)

population <- read.csv(populationCountFile)

# Some data reshaping functions
get_count_df <- function(df_row) {
  if(nrow(df_row) == 1) {
    start_date <- strptime('2020-01-01',format='%Y-%m-%d',tz="GMT")
    n <- names(df_row)
    M <- ncol(df_row) - 4
    count_series <- vector("integer", M)
    date_series <- vector("integer", M)
    for(i in 5:(5+M-1)) {
      count_series[i-4] <- df_row[[i]]
      date_series[i-4] <- difftime(strptime(n[i],format='X%m.%d.%y',tz="GMT"),start_date,units='days')+1
    }
    return(data.frame(count=count_series, dayNum=date_series))
  } else
    return(NULL)
}

# N is the population of the country
get_site_df <- function(province,country, N) {
  conf_df <- get_count_df(jhu_confirmed[jhu_confirmed$Country.Region==country & jhu_confirmed$Province.State==province,])
  if(is.null(conf_df))
    return(NULL);
  names(conf_df)[1] <- 'Confirmed'
  rec_df <- get_count_df(jhu_recovered[jhu_recovered$Country.Region==country & jhu_recovered$Province.State==province,])
  if(is.null(rec_df))
    return(NULL);
  names(rec_df)[1] <- 'Recovered'
  death_df <- get_count_df(jhu_deaths[jhu_deaths$Country.Region==country & jhu_deaths$Province.State==province,])
  if(is.null(death_df))
    return(NULL);
  names(death_df)[1] <- 'Deaths'
  
  
  
  merged1 <- merge(conf_df,rec_df)
  merged2 <- merge(merged1,death_df)
  
  merged2$Removed <- merged2$Deaths + merged2$Recovered
  merged2$Infected <- merged2$Confirmed - merged2$Removed
  merged2$Susceptible <- N - merged2$Infected - merged2$Removed
  #print(merged2)
  return(merged2)
}

globalTs <- NULL
globalPopulation <- 7700000000

for(i in (1:nrow(population))) {
  province <- as.character(population[i,1])
  country <- as.character(population[i,2])
  popCount <- population[i,3]
  #print(paste0("Extracting province '",province,"' state '",country,"' with population ",popCount))
  df <- get_site_df(province,country,popCount)
  if(is.null(df)) {
    print(paste0("Can't find data for province '",province,"' country '",country,"'"))
    next;
  }
  if(is.null(globalTs))
    globalTs <- df
  else
    globalTs[,2:6] <- globalTs[,2:6] + df[,2:6]
  
  maxConfirmed <- max(df$Confirmed)
  maxRemoved <- max(df$Removed)
  if(maxConfirmed < minAllowedConfirmedCount) {
    print(paste0("Skipping province '",province,"' country '",country,"' due to too low confirmed cases: ",maxConfirmed))
    next; 
  }
  if(maxRemoved < minAllowedRemovedCount) {
    print(paste0("Skipping province '",province,"' country '",country,"' due to too low removed cases: ",maxRemoved))
    next; 
  }
  infMax <- max(df$Infected)
  #print(df$Infected/infMax*10.0)
  infCoarse <- round(df$Infected/infMax*30.0)
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
  
  if(nchar(province)>0)
    province <- paste0(province,'@')
  key <- paste0(province, country)
  key <- str_replace_all(key," ","_")
  key <- str_replace_all(key,"\\*","#")
  outFile <- file.path(outDir,paste0(key,'.csv'))
  #print(paste0('Writing ',outFile))
  write.csv(df,outFile, row.names = F)
}

#writing global
globalTs$Susceptible <- globalPopulation - globalTs$Infected - globalTs$Removed
write.csv(globalTs,globalOutputFile, row.names = F)

print("Done")
