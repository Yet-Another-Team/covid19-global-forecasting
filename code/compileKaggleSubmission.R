require(stringr)
require(pracma)

args = commandArgs(trailingOnly=TRUE)

kaggleTestFile <- args[1]
outFile <- args[2]
confirmedCasesSigmoidParamsFile <- args[3]
deathRatesParamsFile <- args[4]

encodeChars <- function(key) {
  key <- str_replace_all(key," ","_")
  key <- str_replace_all(key,"\\*","#")
  key
}

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
kaggleKeys <- data.frame(Province_State=kaggleKeysDf$Province_State, Country_Region=kaggleKeysDf$Country_Region)

remainingKeys <- kaggleKeys

# pre-fitted models predictons
deathRatesParamsDf <- read.csv(deathRatesParamsFile)
confirmedCasesSigmoidParamsDf <- read.csv(confirmedCasesSigmoidParamsFile)

confirmedCasesSigmoidParamsDf$key <- getKey(confirmedCasesSigmoidParamsDf$Province, confirmedCasesSigmoidParamsDf$Country)

getFatalities <- function(dayNums,ConfirmedCases,lag,fatalityRate) {
  maxConfirmed <- max(ConfirmedCases)
  confirmedFunction <- approxfun(dayNums,ConfirmedCases, yleft =0,yright = maxConfirmed)
  
  confirmedFunction(dayNums - lag)*fatalityRate
}

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

ConfirmedSigmoidPredictor <- list(
  name='Confirmed cases sigmoid predictor',
  fun=function(province,country,dates) {
    print(paste0('predicting province ',province,' country ',country))
    province <- encodeChars(province)
    country <- encodeChars(country)
    
    N <- length(dates)
    start_date <- strptime('01-01-2020',format='%m-%d-%Y',tz="GMT")
    dayNums <- difftime(strptime(dates,format='%Y-%m-%d',tz="GMT"),start_date,units='days')+1
    prediction <- data.frame(ConfirmedCases=rep(0,N),Fatalities=rep(0,N))
    
    confSigParamsRow <- confirmedCasesSigmoidParamsDf[(confirmedCasesSigmoidParamsDf$Province == province) & (confirmedCasesSigmoidParamsDf$Country == country),]
    stopifnot(nrow(confSigParamsRow) == 1)
    
    deathRateParamsRow <- deathRatesParamsDf[(confirmedCasesSigmoidParamsDf$Province == province) & (confirmedCasesSigmoidParamsDf$Country == country),]
    stopifnot(nrow(deathRateParamsRow) == 1)
    
    inflDayNum <- confSigParamsRow$ConfirmedInflectionDayNum
    k <- confSigParamsRow$ConfirmedK
    suscPop <- confSigParamsRow$EstimatedSusceptiblePopulation
    
    print(paste0('province ',province,' country ',country,' inflDayNum ',inflDayNum,' k ',k,' suscPop ',suscPop))
    
    deathRateLag <- deathRateParamsRow$realDaysLag
    fatalityRate <- deathRateParamsRow$fatalityRate

    prediction[,1] <- sigmoid(as.numeric((dayNums - inflDayNum)/k))*suscPop
    #print(paste0("confirmed ",prediction[,1]))
    prediction[,2] <- getFatalities(dayNums,prediction$ConfirmedCases,deathRateLag,fatalityRate)
    prediction
  }
)

predictors <- list(
  #jhuTsSIRPredictor,
  ConfirmedSigmoidPredictor)

predictorIdx <- 1
while ((nrow(remainingKeys)>0) && (predictorIdx <= length(predictors)) ) {
  print(paste0('There are ',nrow(remainingKeys),' locations (as defined in Kaggle Challendge) still to predict'))
  predictor <- predictors[[predictorIdx]]
  predictFun <- predictor$fun
  print(paste0("Currently using prediction source '",predictor$name,"'"))
  curPredictedKeys <- c()
  N <- nrow(remainingKeys)
  for(i in (1:N)) {
    province <- as.character(remainingKeys[i,1])
    country <- as.character(remainingKeys[i,2])
    
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