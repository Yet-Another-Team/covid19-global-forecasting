args = commandArgs(trailingOnly=TRUE)

require(stringr)

jhuInputDir <- args[1]
populationCountFile <- args[2]
outputPattern <- args[3] # should end with *.csv

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
  
  if(nchar(province)>0)
    province <- paste0(province,'@')
  key <- paste0(province, country)
  key <- str_replace_all(key," ","_")
  key <- str_replace_all(key,"\\*","#")
  outFile <- file.path(outDir,paste0(key,'.csv'))
  #print(paste0('Writing ',outFile))
  write.csv(df,outFile, row.names = F)
}

print("Done")
