require(stringr)

args = commandArgs(trailingOnly=TRUE)

kaggleTrFile <- args[1]
kaggleTestFile <- args[2]
outFile <- args[3]
jhuTsBasedTrainedSIRdir <- args[4]
lethalRatesFile <- args[5]

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

kaggleTupleFromKey <- function(key) {
  atPos = regexpr('@',key)
  province <- vector(mode='character',length = length(key))
  countryPos <- vector(mode='integer',length = length(key))
  resDf <- data.frame(province=province,countryPos = countryPos)
  resDf$countryPos <- 1
  resDf$countryPos[atPos >0] <- atPos[atPos>0]+1
  resDf$country <- substr(key,resDf$countryPos,nchar(key))
  resDf$country <- str_replace_all(resDf$country,"_"," ")
  resDf$country <- str_replace_all(resDf$country,"#","\\*")
  resDf$province <- as.character(resDf$province)
  resDf$province[atPos>0] <- substr(key[atPos>0],rep(1,sum(atPos>0)),atPos[atPos>0]-1)
  resDf$province <- str_replace_all(resDf$province,"_"," ")
  resDf$province <- str_replace_all(resDf$province,"#","\\*")
  names(resDf) <- c('Province_State','atPos','Country_Region')
  resDf[,c(1,3)]
}

kaggleTestDf <- read.csv(kaggleTestFile)
kaggleTestDf$ConfirmedCases <- 0
kaggleTestDf$Fatalities <- 0

kaggleKeysDf <- kaggleTestDf[,2:3]
kaggleKeysDf <- kaggleKeysDf[!duplicated(kaggleKeysDf),]
kaggleKeys <- data.frame(kaggleKeysDf$Province_State, kaggleKeysDf$Country_Region)

remainingKeys <- kaggleKeys

# pre-fitted SIR models predictons
lethalRatesDf <- read.csv(lethalRatesFile)
lethalRatesDf$Key <- getKey(lethalRatesDf$Province,lethalRatesDf$Country)

predictSIR <- function(predictionsFile,lethalRate,dates) {
  N <- length(dates)
  #print(paste0('opening ',predictionsFile))
  predDf <- read.csv(predictionsFile)
  prediction <- data.frame(Infected=rep(0,N),Removed=rep(0,N))
  start_date <- strptime('01-01-2020',format='%m-%d-%Y',tz="GMT")
  minDay <- min(predDf$days)
  maxDay <- max(predDf$days)
  for(i in (1:N)) {
    curDate <- dates[i]
    cur_date <- strptime(curDate,format='%Y-%m-%d',tz="GMT")
    dayNum <- difftime(cur_date,start_date,units='days')+1
    if(dayNum<minDay)
      prediction[i,] <- predDf[predDf$days == minDay,5:6]
    else  if(dayNum>maxDay) {
      prediction[i,] <- predDf[predDf$days == maxDay,5:6]
    } else {
      prediction[i,] <- predDf[predDf$days == dayNum,5:6]
    }
  }
  prediction$ConfirmedCases <- prediction[,1] + prediction[,2] # infected , removed
  prediction$Fatalities <- round(prediction[,2]*lethalRate)
  prediction$ConfirmedCases <- round(prediction$ConfirmedCases)
  prediction[,3:4]
}

jhuTsSIRPredictor <- list(
  name="JHU timeseries data based fitted SIR models",
  fun=function(province,country,dates) {
    key <- getKey(province = province, country = country)
    predictionFile <- file.path(jhuTsBasedTrainedSIRdir,'per_location_prediction',paste0(key,'.csv'))
    if(!file.exists(predictionFile))
      return(NULL)
    else {
      lr_row <- lethalRatesDf[lethalRatesDf$Key == key,]
      if(nrow(lr_row) == 0) {
        print(paste0("WARNING: can't find lethal rates for key ",key))
        return(NULL)
      }
      lr <- lr_row$DeathToRemovedRatio
      return(predictSIR(predictionFile,lr, dates))
    }
  }
)

# Kaggle Train Set const value extrapolation prediction source
kaggleTrDf <- read.csv(kaggleTrFile)

kaggleTrainSetConstExtrapolationPredictor <- list(
  name='Kaggle Train Set const value extrapolation',
  fun=function(province,country,dates) {
    print(paste0("WARNING: Using train set const-extrapolation as there are no other models for province '",province,"' of country '",country,"' available"))
    N <- length(dates)
    prediction <- data.frame(ConfirmedCases=rep(0,N),Fatalities=rep(0,N))
    curLocTrSet <- kaggleTrDf[(kaggleTrDf$Province_State == province) & (kaggleTrDf$Country_Region == country),]
    start_date <- strptime('01-01-2020',format='%m-%d-%Y',tz="GMT")
    cur_date <- strptime(curLocTrSet$Date,format='%Y-%m-%d',tz="GMT")
    curLocTrSet$DayNum <- difftime(cur_date,start_date,units='days')+1
    minDay <- min(curLocTrSet$DayNum)
    maxDay <- max(curLocTrSet$DayNum)
    for(i in (1:N)) {
      curDate <- dates[i]
      cur_date <- strptime(curDate,format='%Y-%m-%d',tz="GMT")
      dayNum <- difftime(cur_date,start_date,units='days')+1
      if(dayNum<minDay)
        prediction[i,] <- curLocTrSet[curLocTrSet$DayNum == minDay,5:6]
      else  if(dayNum>maxDay) {
        prediction[i,] <- curLocTrSet[curLocTrSet$DayNum == maxDay,5:6]
      } else {
        prediction[i,] <- curLocTrSet[curLocTrSet$DayNum == dayNum,5:6]
      }
    }
    prediction
  }
)

predictors <- list(
  jhuTsSIRPredictor,
  kaggleTrainSetConstExtrapolationPredictor)

predictorIdx <- 1
while ((nrow(remainingKeys)>0) && (predictorIdx <= length(predictors)) ) {
  print(paste0('There are ',nrow(remainingKeys),' locations (as defined in Kaggle Challendge) still to predict'))
  predictor <- predictors[[predictorIdx]]
  predictFun <- predictor$fun
  print(paste0("Currently using prediction source '",predictor$name,"'"))
  curPredictedKeys <- c()
  N <- nrow(remainingKeys)
  for(i in (1:N)) {
    province <- remainingKeys[i,1]
    country <- remainingKeys[i,2]
    
    curPredIndicator <- (kaggleTestDf$Province_State == province) & (kaggleTestDf$Country_Region == country)
    curPredIndices <- which(curPredIndicator)
    testDfAtCurLoc <- kaggleTestDf[curPredIndices,]
    #dates needed to be predicted
    datesToPredict <- as.character(testDfAtCurLoc$Date) # dates in format YYYY-mm-dd
    # trying to predict using current predictor
    predictions <- predictFun(province,country,datesToPredict)
    if(is.null(predictions)) {
      # predictor failed to provide data for current province/country. Fo nothing
      next;
    }
    curPredictedKeys <- c(curPredictedKeys,i)
    kaggleTestDf[curPredIndices,]$ConfirmedCases <- as.integer(predictions$ConfirmedCases)
    kaggleTestDf[curPredIndices,]$Fatalities <- as.integer(predictions$Fatalities)
  }
  print(paste0(length(curPredictedKeys)," locations successfuly predicted using '",predictor$name,"'"))
  remainingKeys <- remainingKeys[-curPredictedKeys,]
  predictorIdx <- predictorIdx + 1
}
res <- kaggleTestDf[,c(1,5,6)]
write.csv(res,file=outFile,row.names = F,quote = F)
print('Done')