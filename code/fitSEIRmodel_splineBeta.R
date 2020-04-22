args = commandArgs(trailingOnly=TRUE)

require(RJSONIO)
require(ggplot2)
require(tidyr)
require(pracma)
require(deSolve)
require(gridExtra)

sirObsFile <- args[1]
pretrainedFile <- args[2]
paramsOutFile <- args[3]
pngOutFile <- args[4]
predTableOutFile <- args[5]

betaChangeIntervalDays <- 7

start_date <- as.Date(strptime('01-01-2020',format='%m-%d-%Y',tz="GMT"))

fineIters <- 2

sigma <- 1.0 / 3.0 # 3 days

softplus.shift <- log(exp(1)-1)
softplus <- function(x) {
  res <- x
  validPos <- (!is.na(x)) & (x<50)
  res[validPos] <- log(1+exp(x[validPos])) # to avoid Inf/
  res
}
# transforms (-inf,inf) -> (0, inf) while equals to 1.0 at 0.0
nonNegativeParam <- function(x) {softplus(x + softplus.shift)}

SEIR <- function(t, state, parameters, get_beta) {
  with(as.list(c(state, parameters)),{
    
    beta <- get_beta(t)
    # rate of change
    
    dS = -beta*I*S/N
    dE = beta*I*S/N - sigma*E
    dI = sigma*E - gamma*I
    dR = gamma*I
    
    # return the rate of change
    list(c(dS, dE, dI, dR))
  })
}

obsSource <- read.csv(sirObsFile)

popCount <- obsSource$Susceptible[1]

maxConfirmed <- max(obsSource$Confirmed)

#print(paste0("popCount ",popCount))
obs <- data.frame(
  days=obsSource$dayNum,
  infected.obs=obsSource$Infected,
  removed.obs=obsSource$Removed)

firstDayIdx <- obsSource$dayNum[which(obsSource$Infected>0)[1]]

#print(paste0('firstDay real Idx ',firstDayIdx))

getPrediction <- function(p) {
  zeroDayNum = p[["zeroDayNum"]]
  firstDayInfectedCount = nonNegativeParam(p[["firstDayInfectedCount"]])
  firstDayExposedCount <- nonNegativeParam(p[["firstDayExposedCount"]])
  #sigma = nonNegativeParam(p[["sigma"]])
  gamma = sigmoid(p[["gamma"]]) # fixed fraction "gamma"  of the infected group will recover during any given day
  N = popCount
  
  b1_idx <- which(names(p) == 'b1')
  paramsCount <- length(p)
  
  # beta spline model construction
  effective_b_params <- unlist(p[b1_idx:paramsCount])
  nonNegBetaScales <- nonNegativeParam(effective_b_params)
  bCount <- length(effective_b_params)
  get_beta <- NULL
  if(bCount>1) {
    betaKnots <- rep(0,bCount) # multipling all of the passed weeks b_scales
    knotDayNum <- rep(0,bCount)
    for(i in 1:bCount){
      proIndicator <- c(rep(T,i),rep(F,(bCount-i)))
      betaKnots[i] <- prod(nonNegBetaScales[proIndicator])
      knotDayNum[i] <- zeroDayNum + betaChangeIntervalDays*i - (betaChangeIntervalDays/2) # middle of the week
    }
    
    # forcing spline to beconstant at the right bound
    betaKnots <- c(betaKnots, rep(betaKnots[bCount],5))
    knotDayNum <- c(knotDayNum, seq(knotDayNum[bCount] + betaChangeIntervalDays, knotDayNum[bCount] + 5*betaChangeIntervalDays, by=betaChangeIntervalDays))
    
    get_beta <- splinefun(knotDayNum, betaKnots, method = "natural")
  } else {
    get_beta <- function(t) nonNegBetaScales
  }
  
  odeDays <- seq(zeroDayNum, zeroDayNum+365, by = 0.1)
  
  
  odeParameters <- c(gamma=gamma, N=N)
  #odeParameters <- c(odeParameters,unlist(p[b1_idx:paramsCount])) # mixing up b1 .. bN
  odeZeroState <- c(S=N,E=firstDayExposedCount, I = firstDayInfectedCount,R=0)
  
  out <- as.data.frame(ode(
    y = odeZeroState,
    times = odeDays,
    func = SEIR,
    parms = odeParameters,
    get_beta = get_beta
    ))
  names(out) <- c('days','susceptible.pred','exposed.pred','infected.pred','removed.pred')
  #print(out)
  
  firstOdeDay <- ceiling(out$days[1])
  lastOdeDay <- floor(out$days[length(out$days)])
  #day aligned
  outAligned <- data.frame(days= firstOdeDay:lastOdeDay)
  
  outAligned$susceptible.pred <- approxfun(out$days,out$susceptible.pred)(outAligned$days)
  outAligned$exposed.pred <- approxfun(out$days,out$exposed.pred)(outAligned$days)
  outAligned$infected.pred <- approxfun(out$days,out$infected.pred)(outAligned$days)
  outAligned$removed.pred <- approxfun(out$days,out$removed.pred)(outAligned$days)
  
  outAligned
}

getPredRunTable<-function(p) {
  out <- getPrediction(p)
  
  outMaxExpV <- max(out$exposed.pred)
  outMaxExpDay <- out$days[which(out$exposed.pred == outMaxExpV)[1]] 
  
  
  outMaxV <- max(out$infected.pred)
  outMaxDay <- out$days[which(out$infected.pred == outMaxV)[1]]
  
  obsCur <- obs
  
  m1 <- merge(obsCur,out, by='days',all.x = T)
  
  # filling up out of range values
  erliestOut <- out[1,]
  latestOut <- out[nrow(out),]
  
  earliestSuscPred <- out$susceptible.pred[1]
  earliestExpPred <- out$exposed.pred[1]
  
  if(sum(is.na(m1$susceptible.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$susceptible.pred) & m1$days<=erliestOut$days,]$susceptible.pred <- erliestOut$susceptible.pred
  if(sum(is.na(m1$susceptible.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$susceptible.pred) & m1$days>=latestOut$days,]$susceptible.pred <- latestOut$susceptible.pred
  
  if(sum(is.na(m1$exposed.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$exposed.pred) & m1$days<=erliestOut$days,]$exposed.pred <- erliestOut$exposed.pred
  if(sum(is.na(m1$exposed.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$exposed.pred) & m1$days>=latestOut$days,]$exposed.pred <- latestOut$exposed.pred
  
  
  if(sum(is.na(m1$infected.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$infected.pred) & m1$days<=erliestOut$days,]$infected.pred <- 0
  if(sum(is.na(m1$infected.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$infected.pred) & m1$days>=latestOut$days,]$infected.pred <- latestOut$infected.pred
  
  if(sum(is.na(m1$removed.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$removed.pred) & m1$days<=erliestOut$days,]$removed.pred <- 0
  if(sum(is.na(m1$removed.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$removed.pred) & m1$days>=latestOut$days,]$removed.pred <- latestOut$removed.pred
  
  
  res <- list(table =m1, peakDay = outMaxDay, peakHeight = outMaxV, peakExpDay= outMaxExpDay, peakExpHeight = outMaxExpV)
}

rmse <- function(obs,pred) {
  d1 <- obs-pred
  sqrt(mean(d1*d1))
}

rmsle <- function(obs,pred) {
  log_obs <- log(max(1,obs+1))
  log_pred <- log(max(1,pred+1))
  sqrt(mean((log_obs-log_pred)*(log_obs-log_pred)))
}

obs_match_loss <- function(loss,
                          infected.obs,
                          removed.obs,
                          infected.pred,
                          removed.pred) {
  (
      loss(infected.obs,infected.pred) * 0.7 +
      loss(removed.obs,removed.pred) * 0.3
  )
}

toMinimize <- function(p) {
  pred <- getPredRunTable(p)
  m1 <- pred$table
  
  l1 <- obs_match_loss(
    rmse,
    m1$infected.obs,
    m1$removed.obs,
    m1$infected.pred,
    m1$removed.pred)
  
  b1_idx <- which(names(p) == 'b1')
  penalty_b_start_idx <- b1_idx + 1
  paramsCount <- length(p)
  betaChangePenalty <- 0
  if(penalty_b_start_idx <= paramsCount) {
    pRange <- unlist(p[penalty_b_start_idx:paramsCount])
    betaChangePenalty <- mean(pRange*pRange) # L2 norm
  }
  
  startAlignmentPenalty <- abs(p[["zeroDayNum"]] - firstDayIdx)
  
  loss <- l1 * (1 + betaChangePenalty*1e-1 + startAlignmentPenalty*1e-2)
  return(loss)
}


# how many weeks to model?
bParamsCount <- ((max(obsSource$dayNum) - min(obsSource$dayNum) + 1) %/% betaChangeIntervalDays)+1 # weeks
print(paste0(bParamsCount," modelling intervals to model"))


pretrainedPar <- fromJSON(pretrainedFile)

#print("Loaded pretrained fit:")
startPar <- list(
  gamma=pretrainedPar[["Gamma"]],
  zeroDayNum = pretrainedPar[["FirstDayNum"]],
  firstDayInfectedCount = pretrainedPar[["FirstDayInfectedCount"]],
  firstDayExposedCount = 0 # this will be reset on each train run setup
)
bRaw <- pretrainedPar[["bRaw"]]
if(length(bRaw) < bParamsCount)
  bRaw <- c(bRaw, rep(0, bParamsCount - length(bRaw)))
names(bRaw) <- paste0("b",(1:(length(bRaw))))
startPar <- c(startPar,bRaw)
#print(startPar)


print("Final precise fit began")
if(bParamsCount == 1) {
  fineIters <- 2
}

optFineRes <- NULL
firstBparamIdx <- which(names(startPar)=="b1")
totalParamCount <- length(startPar)
# long train run
i <- 1
while(i <= fineIters) {
  set.seed(12323 + 123*i)
  bDeltas <- rnorm(bParamsCount,mean=0,sd = 1e-3 * max(1,i-10))
  fineStartPars <- startPar
  fineStartPars$firstDayExposedCount <- runif(1,min = 0,max= 10)
  fineStartPars[firstBparamIdx:totalParamCount] <-
    unlist(fineStartPars[firstBparamIdx:totalParamCount]) + bDeltas
  
  methodToUse <- switch((i %% 3) + 1,
                        "CG",
                        "BFGS",
                        "Nelder-Mead"
                        )
  derivativeBasedMaxit <- as.integer(900+bParamsCount*100)
  nmMaxit <- as.integer(9000+bParamsCount*1000)
  maxit <- switch((i %% 3) +1,
                  derivativeBasedMaxit,
                  derivativeBasedMaxit,
                  nmMaxit
                  )
  
  optCtr <- list(trace=0,maxit=maxit) # set trace to value higher than 0, if you want details
  
  options(warn=-1)
  sink("NUL")
  
  #finRes <- optim(fineStartPars, toMinimize,control = optCtr, method=methodToUse)
  
  finRes <- tryCatch(
    optim(fineStartPars, toMinimize,control = optCtr, method=methodToUse),
    error = function(e) NULL
      )
  
  sink()
  options(warn=0)
  
  if(!is.null(finRes) && (is.null(optFineRes) || optFineRes$value > finRes$value)) {
    print(paste0("Iteration ",i," (",methodToUse,"): ----> fine loss improved from ",optFineRes$value," to ",finRes$value))
    optFineRes <- finRes
  } else {
    print(paste0("Iteration ",i," (",methodToUse,"): ",finRes$value," (best is still ",optFineRes$value,")"))
  }
  if(is.null(finRes)) {
    fineIters <- fineIters+1
  }
  
  i <- i + 1
}

# coersing beta scale for future to 0.0 (thus scales to 1.0)
b1_idx <- which(names(optFineRes$par) == 'b1')
paramsCount <- length(optFineRes$par)

zeroDayNum <- optFineRes$par[["zeroDayNum"]]
latestObsDayNum <- max(obsSource$dayNum)
weeksIncluded <- as.integer(ceiling((latestObsDayNum - zeroDayNum + 1) / betaChangeIntervalDays))
firstParToZero <- b1_idx + weeksIncluded
if(firstParToZero <= paramsCount) {
  print(paste0("Coersing ",(paramsCount-firstParToZero+1)," future beta scales"))
  optFineRes$par <- optFineRes$par[1:(firstParToZero-1)]
}
bestPrediction <- getPredRunTable(optFineRes$par)


# now calculating beta

bScales <- nonNegativeParam(unlist(optFineRes$par[b1_idx:(length(optFineRes$par))]))
names(bScales) <- NULL
bCount <- length(bScales)
beta <- rep(0.0,bCount)
for(i in 1:bCount){
  proIndicator <- c(rep(T,i),rep(F,(bCount-i)))
  beta[i] <- prod(bScales[proIndicator])
}

print("Fine optimization finished. Results")
#print(optFineRes)

rmse = optFineRes$value
zeroDayNum = optFineRes$par[["zeroDayNum"]]
firstDayInfectedCount = nonNegativeParam(optFineRes$par[["firstDayInfectedCount"]])
#sigma = nonNegativeParam(optFineRes$par[["sigma"]])
gamma = sigmoid(optFineRes$par[["gamma"]])
r0 = beta/gamma

recentBeta <- beta[length(beta)]
recentR0<- r0[length(beta)]

paramsList <- list()
paramsList$Gamma <- gamma
paramsList$Sigma <- sigma
paramsList$WeeklyBeta <- beta
paramsList$RecentBeta <- recentBeta
if(length(bScales)>=2) {
  paramsList$BetaScales <- bScales[2:length(bScales)]
}
paramsList$WeeklyR0 <- r0
paramsList$RecentR0 <- recentR0
paramsList$FirstDayNum <- zeroDayNum
paramsList$FirstDayInfectedCount <- firstDayInfectedCount
paramsList$PeakDayNum <- bestPrediction$peakDay
paramsList$PeakDayInfectedCount <- bestPrediction$peakHeight
paramsList$TotalPopulation <- popCount
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
    (predTableG$group == 'infected.pred') |
      (predTableG$group == 'removed.pred') |
      (predTableG$group == 'susceptible.pred') |
      (predTableG$group == 'exposed.pred')  ,]$Type <- 'Model'
  predTableG$Type <- as.factor(predTableG$Type)
  
  # adding Group feature: susceptible / infected / removed
  predTableG$Group <- 'Susceptible'
  predTableG[(predTableG$group == 'infected.obs') | (predTableG$group == 'infected.pred'),]$Group <- 'Infected'
  predTableG[(predTableG$group == 'removed.obs') | (predTableG$group == 'removed.pred'),]$Group <- 'Removed'
  predTableG[(predTableG$group == 'exposed.pred'),]$Group <- 'Exposed'
  predTableG$Group <- as.factor(predTableG$Group)
  
  if(sum(predTableG$people == 0)>0)
    predTableG[predTableG$people == 0,]$people <- 1e-5 # to avoid problems with log scale
  
  predTableG <- predTableG[predTableG$Group != 'Susceptible',]
  
  obsOnlyG <- predTableG[predTableG$Type == 'Actual',]
  
  cols <- c("Infected" = "red","Exposed"="orange", "Removed" = "green", "Susceptible"="blue")
  shapes <- c("factor"=1)
  
  obsOnlyG$Date <- start_date + obsOnlyG$days - 1
  
  res <-
    p +
    geom_point(aes(x=Date,y=people,fill=Group),data=obsOnlyG, shape=21,color='transparent') + 
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
    (predTableG$group == 'infected.pred') |
      (predTableG$group == 'removed.pred') |
      (predTableG$group == 'susceptible.pred') |
      (predTableG$group == 'exposed.pred')  ,]$Type <- 'Model'
  predTableG$Type <- as.factor(predTableG$Type)
  
  # adding Group feature: susceptible / infected / removed
  predTableG$Group <- 'Susceptible'
  predTableG[(predTableG$group == 'infected.obs') | (predTableG$group == 'infected.pred'),]$Group <- 'Infected'
  predTableG[(predTableG$group == 'removed.obs') | (predTableG$group == 'removed.pred'),]$Group <- 'Removed'
  predTableG[(predTableG$group == 'exposed.pred'),]$Group <- 'Exposed'
  predTableG$Group <- as.factor(predTableG$Group)
  
  modelOnlyG <- predTableG[predTableG$Type == 'Model',]

  cols <- c("Infected" = "red","Exposed"="orange", "Removed" = "green", "Susceptible"="blue")

  modelOnlyG$Date <- start_date + modelOnlyG$days - 1
  
  res <- p +
    scale_color_manual(values = cols,name="Model Prediction") +
    geom_line(aes(x=Date,y=people,color=Group),data=modelOnlyG,size=0.5)
    
  #theme(legend.position = "bottom")
  res
}


max_val <- max(c(obsSource$Infected,obsSource$Removed)) * 1.1
latest_obs <- max(obsSource$dayNum) + 1
earliest_obs <- min(obsSource$dayNum)

obsFName <- basename(sirObsFile)
obsFName <- substr(obsFName,1,(nchar(obsFName)-4))
atPos = regexpr('@',obsFName)[1]
if(popCount > 5e9) {
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
  labs(title = paste0("SEIR model fit [",descr,"]"),
       subtitle = paste0("recent R0 = ",round(recentR0,1),
              " recent beta = ",round(recentBeta,3),
              ' sigma = ',round(sigma,3),
              ' gamma = ',round(gamma,3),
              ' loss = ',round(rmse)))+
  scale_y_continuous(limits=c(0,max_val))
  #scale_x_continuous(limits=c(earliest_obs,latest_obs))

p2 <- ggplot()
yearPred <- getPrediction(optFineRes$par)
max_pred_val <- max(c(obsSource$Infected,obsSource$Removed,yearPred$infected.pred,yearPred$removed.pred, yearPred$exposed.pred)) * 1.1
p2 <- plotObsTable(bestPrediction$table,p2)
p2 <- plotPredTable(yearPred,p2)
p2 <- p2 + labs(title = "One year simulation",
                caption = "COVID-19 epidemic dynamics model") +
  guides(fill=FALSE, color=FALSE) +
  scale_y_continuous(limits=c(0,max_pred_val))


p3 <- grid.arrange(p, p2, nrow=2)

ggsave(pngOutFile,p3)
print("Figure saved")

predCols <- ncol(yearPred)
yearPred$Date <- as.character(as.Date(strptime('2020-01-01',format='%Y-%m-%d',tz="GMT"))+yearPred$days-1)
yearPred <- yearPred[,c(predCols+1,1:predCols)]
write.csv(yearPred, file=predTableOutFile, row.names = F)
print("Predicion table written")
print("Done")