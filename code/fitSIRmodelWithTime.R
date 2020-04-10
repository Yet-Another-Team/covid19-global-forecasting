args = commandArgs(trailingOnly=TRUE)

require(RJSONIO)
require(ggplot2)
require(tidyr)
require(pracma)
require(deSolve)
require(gridExtra)

sirObsFile <- args[1]
paramsOutFile <- args[2]
pngOutFile <- args[3]
predTableOutFile <- args[4]

betaChangeIntervalDays <- 7

optIters <- 40
fineIters <- 5

softplus.shift <- log(exp(1)-1)
softplus <- function(x) {
  res <- x
  validPos <- (!is.na(x)) & (x<50)
  res[validPos] <- log(1+exp(x[validPos])) # to avoid Inf/
  res
}
# transforms (-inf,inf) -> (0, inf) while equals to 1.0 at 0.0
nonNegativeParam <- function(x) {softplus(x + softplus.shift)}


SIR <- function(t, state, parameters, get_beta) {
  with(as.list(c(state, parameters)),{
    
    beta <- get_beta(t,parameters)
    # rate of change
    dS = - beta *I*S/N
    dI = beta*I*S/N - gamma*I
    dR = gamma*I
    
    # return the rate of change
    list(c(dS, dI, dR))
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
  gamma = sigmoid(p[["gamma"]]) # fixed fraction "gamma"  of the infected group will recover during any given day
  N = popCount
  zeroDayNum <- round(zeroDayNum,1) # to align with 0.1 simulation step
  
  b1_idx <- which(names(p) == 'b1')
  paramsCount <- length(p)
  
  get_beta <- function(t,p) {
    t0 <- t - zeroDayNum
    ode_b1_idx <- 3
    weekIdx <- (t0 %/% betaChangeIntervalDays) # zero based
    weekIdx <- min(weekIdx, (paramsCount - b1_idx))
    
    
    
    # each of the b_params is scale factor K so that b(week_i+1) = b(week_i) * all_b_params(week_num)
    effective_b_params <- unlist(p[ode_b1_idx:(ode_b1_idx+weekIdx)])
    nonNegBetaScales <- nonNegativeParam(effective_b_params)
    beta <- prod(nonNegBetaScales) # multipling all of the passed weeks b_scales
    beta
  }
  
  odeDays <- seq(zeroDayNum, zeroDayNum+365, by = 0.1)
  odeParameters <- c(gamma=gamma, N=N)
  odeParameters <- c(odeParameters,unlist(p[b1_idx:paramsCount])) # mixing up b1 .. bN
  odeZeroState <- c(S=N, I = firstDayInfectedCount,R=0)
  
  out <- as.data.frame(ode(
    y = odeZeroState,
    times = odeDays,
    func = SIR,
    parms = odeParameters,
    get_beta = get_beta
    ))
  names(out) <- c('days','susceptible.pred','infected.pred','removed.pred')
  #print(out)
  
  out <- out[out$days %% 1 ==0,]
  out
}

getPredRunTable<-function(p) {
  out <- getPrediction(p)
  
  outMaxV <- max(out$infected.pred)
  outMaxDay <- out$days[which(out$infected.pred == outMaxV)[1]]
  
  obsCur <- obs
  
  m1 <- merge(obsCur,out, by='days',all.x = T)
  
  # filling up out of range values
  erliestOut <- out[1,]
  latestOut <- out[nrow(out),]
  
  earliestSuscPrep <- out$susceptible.pred[1]
  
  if(sum(is.na(m1$susceptible.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$susceptible.pred) & m1$days<=erliestOut$days,]$susceptible.pred <- earliestSuscPrep
  if(sum(is.na(m1$susceptible.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$susceptible.pred) & m1$days>=latestOut$days,]$susceptible.pred <- latestOut$susceptible.pred
  
  if(sum(is.na(m1$infected.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$infected.pred) & m1$days<=erliestOut$days,]$infected.pred <- 0
  if(sum(is.na(m1$infected.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$infected.pred) & m1$days>=latestOut$days,]$infected.pred <- latestOut$infected.pred
  
  if(sum(is.na(m1$removed.pred) & m1$days<=erliestOut$days)>0)
    m1[is.na(m1$removed.pred) & m1$days<=erliestOut$days,]$removed.pred <- 0
  if(sum(is.na(m1$removed.pred) & m1$days>=latestOut$days)>0)
    m1[is.na(m1$removed.pred) & m1$days>=latestOut$days,]$removed.pred <- latestOut$removed.pred
  
  
  res <- list(table =m1, peakDay = outMaxDay, peakHeight = outMaxV)
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
  
  # composing penalty for beta changes for the future time (out of observed days)
  maxDayNum <- max(m1$days)
  zeroDayNum = p[["zeroDayNum"]]
  
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
    betaChangePenalty <- sum(pRange*pRange) # L2 norm
  }
  
  
  #print(paste0('est pop ',sigmoid(p[5])*popCount,' max confiremed ', maxConfirmed))
  
  loss <- l1 * (1+ betaChangePenalty/10.0)
  return(loss)
}


optRes <- NULL
# pre train
for(i in (1:optIters)) {
  # will try to fit several times with different seeds
  set.seed(123 + 101*i)
  startP = c(zeroDayNum = firstDayIdx, # when the infection started
             firstDayInfectedCount = runif(1,min=1,max=10), # how many infected on the first day
             gamma=runif(1),
             b1 = rnorm(1)
            )
  #startP <- c(startP,bParams)
  #print(startP)
  #print(toMinimize(startP)) # loss value at start
  optCtr <- list(trace=0,maxit=300) # set trace to value higher than 0, if you want details
  #print(toMinimize(startP))
  curOptRes <- optim(startP, toMinimize,control = optCtr)
  if(curOptRes$convergence > 1)
    next; # we analyze only converged results
  
  if(is.null(optRes) || optRes$value > curOptRes$value) {
    print(paste0("Iteration ",i,": loss improved from ",optRes$value," to ",curOptRes$value))
	  optRes <- curOptRes
  } else {
    #print(paste0("Iteration ",i,": loss did not improve (",curOptRes$value,")"))
  }
}

print("Best sparse fit:")
print(optRes)
print(paste0("Beta: ",nonNegativeParam(optRes[["b1"]])))
print(paste0("Gamma: ",sigmoid(optRes[["gamma"]])))
# how many weeks to model?
bParamsCount <- ((max(obsSource$dayNum) - min(obsSource$dayNum) + 1) %/% betaChangeIntervalDays)+1 # weeks
print(paste0(bParamsCount," modelling intervals to model"))


print("Final precise fit began")
preOptimized <- optRes$par
if(bParamsCount == 1) {
  fineIters <- 1
}

optFineRes <- NULL
# long train run
for(i in 1:fineIters) {
  set.seed(12323 + 123*i)
  bScales <- rnorm(bParamsCount-1,mean=0,sd = 1e-2)
  names(bScales) <- paste0("b",2:bParamsCount)
  fineStartPars <- c(preOptimized,bScales)
  
  optCtr <- list(trace=0,maxit=as.integer(9000+bParamsCount*1000)) # set trace to value higher than 0, if you want details
  #optCtr <- list(trace=1,maxit=1000) # set trace to value higher than 0, if you want details
  finRes <- optim(fineStartPars, toMinimize,control = optCtr)
  if(is.null(optFineRes) || optFineRes$value > finRes$value) {
    print(paste0("Iteration ",i,": loss improved from ",optFineRes$value," to ",finRes$value))
    optFineRes <- finRes
  } else {
    print(paste0("Iteration ",i,": loss did not improve (",finRes$value,")"))
  }
}

bestPrediction <- getPredRunTable(finRes$par)

b1_idx <- which(names(finRes$par) == 'b1')
bScales <- nonNegativeParam(unlist(finRes$par[b1_idx:(length(finRes$par))]))
names(bScales) <- NULL
bCount <- length(bScales)
beta <- rep(0.0,bCount)
for(i in 1:bCount){
  proIndicator <- c(rep(T,i),rep(F,(bCount-i)))
  beta[i] <- prod(bScales[proIndicator])
}

print("Optimization finished. Results")
print(finRes)

rmse = finRes$value
zeroDayNum = finRes$par[["zeroDayNum"]]
firstDayInfectedCount = nonNegativeParam(finRes$par[["firstDayInfectedCount"]])
gamma = sigmoid(finRes$par[["gamma"]])
r0 = beta/gamma

paramsList <- list()
paramsList$Gamma <- gamma
paramsList$Beta <- beta
if(length(bScales)>=2) {
  paramsList$BetaScales <- bScales[2:length(bScales)]
}
paramsList$R0 <- r0
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
  predTableG[(predTableG$group == 'infected.pred') | (predTableG$group == 'removed.pred') | (predTableG$group == 'susceptible.pred'),]$Type <- 'Model'
  predTableG$Type <- as.factor(predTableG$Type)
  
  # adding Group feature: susceptible / infected / removed
  predTableG$Group <- 'Susceptible'
  predTableG[(predTableG$group == 'infected.obs') | (predTableG$group == 'infected.pred'),]$Group <- 'Infected'
  predTableG[(predTableG$group == 'removed.obs') | (predTableG$group == 'removed.pred'),]$Group <- 'Removed'
  predTableG$Group <- as.factor(predTableG$Group)
  
  if(sum(predTableG$people == 0)>0)
    predTableG[predTableG$people == 0,]$people <- 1e-5 # to avoid problems with log scale
  
  predTableG <- predTableG[predTableG$Group != 'Susceptible',]
  
  obsOnlyG <- predTableG[predTableG$Type == 'Actual',]
  
  cols <- c("Infected" = "red", "Removed" = "green", "Susceptible"="blue")
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
  predTableG[(predTableG$group == 'infected.pred') | (predTableG$group == 'removed.pred') | (predTableG$group == 'susceptible.pred'),]$Type <- 'Model'
  predTableG$Type <- as.factor(predTableG$Type)
  
  # adding Group feature: susceptible / infected / removed
  predTableG$Group <- 'Susceptible'
  predTableG[(predTableG$group == 'infected.obs') | (predTableG$group == 'infected.pred'),]$Group <- 'Infected'
  predTableG[(predTableG$group == 'removed.obs') | (predTableG$group == 'removed.pred'),]$Group <- 'Removed'
  predTableG$Group <- as.factor(predTableG$Group)
  
  modelOnlyG <- predTableG[predTableG$Type == 'Model',]

  cols <- c("Infected" = "red", "Removed" = "green", "Susceptible"="blue")

  res <- p +
    scale_color_manual(values = cols,name="Model Prediction") +
    geom_line(aes(x=days,y=people,color=Group),data=modelOnlyG,size=0.5)
    
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
  labs(title = paste0("SIR model fit [",descr,"]"),
       subtitle = paste0("avg R0 = ",round(mean(r0),1),
              " avg beta = ",round(mean(beta),2),
              ' gamma = ',round(gamma,2),
              ' loss = ',round(rmse),
              ' population = ',round(popCount)))+
  scale_y_continuous(limits=c(0,max_val)) +
  scale_x_continuous(limits=c(earliest_obs,latest_obs))

p2 <- ggplot()
yearPred <- getPrediction(optRes$par)
p2 <- plotObsTable(bestPrediction$table,p2)
p2 <- plotPredTable(yearPred,p2)
p2 <- p2 + labs(title = "One year simulation",
                caption = "COVID-19 epidemic dynamics model") +
  guides(fill=FALSE, color=FALSE) +
  scale_y_continuous(limits=c(0,max_val))


p3 <- grid.arrange(p, p2, nrow=2)

ggsave(pngOutFile,p3)
print("Figure saved")

predCols <- ncol(yearPred)
yearPred$Date <- as.character(as.Date(strptime('2020-01-01',format='%Y-%m-%d',tz="GMT"))+yearPred$days-1)
yearPred <- yearPred[,c(predCols+1,1:predCols)]
write.csv(yearPred, file=predTableOutFile)
print("Predicion table written")
print("Done")