args = commandArgs(trailingOnly=TRUE)

require(RJSONIO)
require(ggplot2)
require(tidyr)
require(pracma)
require(gridExtra)

sirObsFile <- args[1]
paramsOutFile <- args[2]
pngOutFile <- args[3]
predTableOutFile <- args[4]

optIters <- 100

obsSource <- read.csv(sirObsFile)

popCount <- obsSource$Susceptible[1]

#print(paste0("popCount ",popCount))
obs <- data.frame(
  days=obsSource$dayNum,
  confirmed.obs=obsSource$Confirmed
  #deaths.obs = obsSource$Deaths
  )

firstDayIdx <- obsSource$dayNum[which(obsSource$Confirmed>0)[1]]

#print(paste0('firstDay real Idx ',firstDayIdx))

getPrediction <- function(p) {
  inflectionDayNum = p[["inflectionDayNum"]]
  k = abs(p[["k"]])+1e-5
  
  #deathInflectionDayNum = p[["deathInflectionDayNum"]]
  #deathK = abs(p[["deathK"]])+1e-5
  
  N = sigmoid(p[["popFactor"]])*popCount
  
  #deathN = sigmoid(p[["deathPopFactor"]])*popCount 
  
  inflectionDayNum <- round(inflectionDayNum)
  
  days <- 1:365
  
  confirmed <- sigmoid((days-inflectionDayNum)/k)*N
  
  #deaths <- sigmoid((days-deathInflectionDayNum)/deathK)*deathN
  
  out <- data.frame(days= days, confirmed.pred = confirmed
                    #, deaths.pred = deaths
                    )
  out
}

getPredRunTable<-function(p) {
  out <- getPrediction(p)
  
  obsCur <- obs
  
  m1 <- merge(obsCur,out, by='days',all.x = T)
  
  latestPredConfirmed = out$confirmed.pred[nrow(out)]
  latestPredDeaths = out$deaths.pred[nrow(out)]
  
  # filling up out of range values
  erliestOut <- out[1,]
  latestOut <- out[nrow(out),]
  
  earliestSuscPrep <- out$susceptible.pred[1]
  
  if(sum(is.na(m1$confirmed.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$confirmed.pred) & m1$days<=erliestOut$days,]$confirmed.pred <- 0
  if(sum(is.na(m1$confirmed.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$confirmed.pred) & m1$days>=latestOut$days,]$confirmed.pred <- latestPredConfirmed
  
  #if(sum(is.na(m1$deaths.pred) & m1$days<=erliestOut$days)>0)
  #  m1[is.na(m1$deaths.pred) & m1$days<=erliestOut$days,]$deaths.pred <- 0
  #if(sum(is.na(m1$deaths.pred) & m1$days>=latestOut$days)>0)
  #  m1[is.na(m1$deaths.pred) & m1$days>=latestOut$days,]$deaths.pred <- latestPredDeaths

  res <- list(table =m1)
}

rmse <- function(obs,pred) {
  sqrt(mean((obs-pred)*(obs-pred)))
}


toMinimize <- function(p) {
  pred <- getPredRunTable(p)
  
  m1 <- pred$table
  
  loss <- rmse(m1$confirmed.pred,m1$confirmed.obs)
    #+ rmse(m1$deaths.pred,m1$deaths.obs)
  return(loss)
}

optRes <- NULL
for(i in (1:optIters)) {
  # will try to fit several times with different seeds
  set.seed(12433 + 101*i)
  startP = c(inflectionDayNum=runif(1)*365, # when the infection started
             k = runif(1)+0.5,
             popFactor = runif(1)#,  # pop factor,
             #deathInflectionDayNum=runif(1)*365, # when the infection started
             #deathK = runif(1)+0.5,
             #deathPopFactor = runif(1)  # pop factor
             ) 
  #print(startP)
  #print(toMinimize(startP)) # loss value at start
  optCtr <- list(trace=0,maxit=10000) # set trace to value higher than 0, if you want details
  curOptRes <- optim(startP, toMinimize,control = optCtr)
  if(curOptRes$convergence != 0)
    next; # we analyze only converged results
  
  if(is.null(optRes) || optRes$value > curOptRes$value) {
    print(paste0("Iteration ",i,": loss improved from ",optRes$value," to ",curOptRes$value))
	optRes <- curOptRes
  } else {
    #print(paste0("Iteration ",i,": loss did not improve (",curOptRes$value,")"))
  }
}

bestPrediction <- getPredRunTable(optRes$par)

rmse = optRes$value
inflectionDayNum = optRes$par[["inflectionDayNum"]]
popFactor = sigmoid(optRes$par[["popFactor"]])
k = abs(optRes$par[["k"]])+1e-5

#deathInflectionDayNum = optRes$par[["deathInflectionDayNum"]]
#deathPopFactor = sigmoid(optRes$par[["deathPopFactor"]])
#deathK = abs(optRes$par[["deathK"]])+1e-5

paramsList <- list()
paramsList$TotalPopulation <- popCount
paramsList$ConfirmedInflectionDayNum <- inflectionDayNum
paramsList$ConfirmedK <- k
paramsList$ConfirmedPopFactor <- popFactor
paramsList$EstimatedSusceptiblePopulation <- popCount*popFactor
#paramsList$TotalDeaths <- popCount*deathPopFactor
#paramsList$DeathsInflectionDayNum <- deathInflectionDayNum
#paramsList$DeathsK <- deathK
#paramsList$DeathsPopFactor <- deathPopFactor

paramsList$Loss <- rmse

exportJson <- toJSON(paramsList)
write(exportJson, paramsOutFile)
print("written param file")

plotObsTable <- function(predTable,p) {
  predTableG <- gather(predTable,key="group",value='people',-days)
  predTableG$group <- as.factor(predTableG$group)
  
  # adding Type feaure : Actual or Model
  predTableG$Type <- 'Actual'
  predTableG[
    (predTableG$group == 'confirmed.pred' |(predTableG$group == 'deaths.pred') ),]$Type <- 'Model'
  predTableG$Type <- as.factor(predTableG$Type)
  
  # adding Group feature: susceptible / infected / removed
  predTableG$Group <- 'Confirmed'
  #predTableG[(predTableG$group == 'deaths.obs') | (predTableG$group == 'deaths.pred'),]$Group <- 'Deaths'
  predTableG$Group <- as.factor(predTableG$Group)
  
  if(sum(predTableG$people == 0)>0)
    predTableG[predTableG$people == 0,]$people <- 1e-5 # to avoid problems with log scale
  
  obsOnlyG <- predTableG[predTableG$Type == 'Actual',]
  
  cols <- c("Confirmed" = "red","Deaths"="grey")
  shapes <- c("factor"=1)
  
  res <-
    p +
    xlab("Days since year 2020 start") +
    geom_point(aes(x=days,y=people,fill=Group),data=obsOnlyG, shape=21,color='transparent') + 
    scale_fill_manual(values = cols ,name="Observations") +
    theme_bw() #+
  res
}

plotPredTable <- function(predTable,p) {
  predTableG <- gather(predTable,key="group",value='people',-days)
  predTableG$group <- as.factor(predTableG$group)
  
  # adding Type feaure : Actual or Model
  predTableG$Type <- 'Actual'
  predTableG[
    (predTableG$group == 'confirmed.pred') |(predTableG$group == 'deaths.pred')  ,]$Type <- 'Model'
  predTableG$Type <- as.factor(predTableG$Type)
  
  # adding Group feature: susceptible / infected / removed
  predTableG$Group <- 'Confirmed'
  #predTableG[(predTableG$group == 'deaths.obs') | (predTableG$group == 'deaths.pred'),]$Group <- 'Deaths'
  predTableG$Group <- as.factor(predTableG$Group)
  
  modelOnlyG <- predTableG[predTableG$Type == 'Model',]

  cols <- c("Confirmed" = "red","Deaths"="grey")

  res <- p +
    scale_color_manual(values = cols,name="Model Prediction") +
    geom_line(aes(x=days,y=people,color=Group),data=modelOnlyG,size=0.5)
    
  #theme(legend.position = "bottom")
  res
}


max_val <- max(obsSource$Confirmed) * 1.1
latest_obs <- max(obsSource$dayNum) + 1
earliest_obs <- min(obsSource$dayNum)

obsFName <- basename(sirObsFile)
obsFName <- substr(obsFName,1,(nchar(obsFName)-4))
atPos = regexpr('@',obsFName)[1]
if(popCount >  6e9) {
  descr <- "worldwide"
} else {
  if(atPos>1)
    province <- paste0(substr(obsFName,1,atPos-1),' - ')
  else
    province <- ''
  country <- substr(obsFName,atPos+1,nchar(obsFName))
  descr <-paste0(province, country)
}

p <- ggplot()
p <- plotObsTable(bestPrediction$table,p)
p <- plotPredTable(bestPrediction$table,p) +
  labs(title = paste0("Sigmoid model fit [",descr,"]"),
       subtitle = paste0("inflection day = ",round(inflectionDayNum),
              " k = ",round(k,2),
              ' loss = ',round(rmse),
              ' init suscept pop = ',round(popFactor*popCount)))+
  scale_y_continuous(limits=c(0,max_val)) +
  scale_x_continuous(limits=c(earliest_obs,latest_obs))

p2 <- ggplot()
yearPred <- getPrediction(optRes$par)
p2 <- plotObsTable(bestPrediction$table,p2)
p2 <- plotPredTable(yearPred,p2)
p2 <- p2 + labs(title = "One year simulation",
                caption = "COVID-19 epidemic dynamics model") +
  guides(fill=FALSE, color=FALSE)


p3 <- grid.arrange(p, p2, nrow=2)

ggsave(pngOutFile,p3)
print("Figure saved")

predCols <- ncol(yearPred)
yearPred$Date <- as.character(as.Date(strptime('2020-01-01',format='%Y-%m-%d',tz="GMT"))+yearPred$days-1)
yearPred <- yearPred[,c(predCols+1,1:predCols)]
write.csv(yearPred, file=predTableOutFile)
print("Predicion table written")
print("Done")