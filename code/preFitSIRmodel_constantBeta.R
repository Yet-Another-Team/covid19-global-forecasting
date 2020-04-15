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


start_date <- as.Date(strptime('01-01-2020',format='%m-%d-%Y',tz="GMT"))

optIters <- 50
iterEvalCount <- 500


softplus.shift <- log(exp(1)-1)
softplus <- function(x) {
  res <- x
  validPos <- (!is.na(x)) & (x<50)
  res[validPos] <- log(1+exp(x[validPos])) # to avoid Inf/
  res
}
# transforms (-inf,inf) -> (0, inf) while equals to 1.0 at 0.0
nonNegativeParam <- function(x) {softplus(x + softplus.shift)}


SIR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
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
  beta = nonNegativeParam(p[["beta"]])
  N = popCount
  
  zeroDayNum <- round(zeroDayNum,1) # to align with 0.1 simulation step
  odeDays <- seq(zeroDayNum, zeroDayNum+365, by = 0.1)
  
  
  odeParameters <- c(gamma=gamma,beta=beta, N=N)
  #odeParameters <- c(odeParameters,unlist(p[b1_idx:paramsCount])) # mixing up b1 .. bN
  odeZeroState <- c(S=N, I = firstDayInfectedCount,R=0)
  
  out <- as.data.frame(ode(
    y = odeZeroState,
    times = odeDays,
    func = SIR,
    parms = odeParameters
    ))
  names(out) <- c('days','susceptible.pred','infected.pred','removed.pred')
  #print(out)
  
  out$days <- round(out$days,1) # as days can disalign with odeDays
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
  
  l1 <- obs_match_loss(
    rmse,
    m1$infected.obs,
    m1$removed.obs,
    m1$infected.pred,
    m1$removed.pred)
  
  lateStartDayPenalty <- max(0, p[["zeroDayNum"]] - firstDayIdx)
  
  l1 + lateStartDayPenalty
}


optRes <- NULL
# pre train
for(i in (1:optIters)) {
  # will try to fit several times with different seeds
  set.seed(123 + 101*i)
  startP = c(zeroDayNum = firstDayIdx, # when the infection started
             firstDayInfectedCount = runif(1,min=1,max=10), # how many infected on the first day
             gamma=runif(1),
             beta = rnorm(1)
            )
  optCtr <- list(trace=0,maxit=iterEvalCount) # set trace to value higher than 0, if you want details
  
  curOptRes <- optim(startP, toMinimize,control = optCtr)
  if(curOptRes$convergence > 1)
    next; # we analyze only converged results
  
  if(is.null(optRes) || optRes$value > curOptRes$value) {
    print(paste0("Iteration ",i,": coarse loss improved from ",optRes$value," to ",curOptRes$value))
	  optRes <- curOptRes
  } else {
    #print(paste0("Iteration ",i,": loss did not improve (",curOptRes$value,")"))
  }
}

print("Best coarse fit:")
#print(optRes)
print(paste0("First day: ",optRes$par[["zeroDayNum"]]))
print(paste0("Beta: ",nonNegativeParam(optRes$par[["beta"]])))
print(paste0("Gamma: ",sigmoid(optRes$par[["gamma"]])))
print(paste0("Loss: ",optRes$value))

rmse = optRes$value
zeroDayNum = optRes$par[["zeroDayNum"]]
firstDayInfectedCount = nonNegativeParam(optRes$par[["firstDayInfectedCount"]])
gamma = sigmoid(optRes$par[["gamma"]])
beta <- nonNegativeParam(optRes$par[["beta"]])
r0 = beta/gamma


exportJson <- toJSON(optRes$par)
write(exportJson, paramsOutFile)
print("written param file")

bestPrediction <- getPredRunTable(optRes$par)


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
  predTableG[(predTableG$group == 'infected.pred') | (predTableG$group == 'removed.pred') | (predTableG$group == 'susceptible.pred'),]$Type <- 'Model'
  predTableG$Type <- as.factor(predTableG$Type)
  
  # adding Group feature: susceptible / infected / removed
  predTableG$Group <- 'Susceptible'
  predTableG[(predTableG$group == 'infected.obs') | (predTableG$group == 'infected.pred'),]$Group <- 'Infected'
  predTableG[(predTableG$group == 'removed.obs') | (predTableG$group == 'removed.pred'),]$Group <- 'Removed'
  predTableG$Group <- as.factor(predTableG$Group)
  
  modelOnlyG <- predTableG[predTableG$Type == 'Model',]

  cols <- c("Infected" = "red", "Removed" = "green", "Susceptible"="blue")

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
  labs(title = paste0("SIR model fit [",descr,"]"),
       subtitle = paste0("R0 = ",round(r0,1),
              " beta = ",round(beta,2),
              ' gamma = ',round(gamma,2),
              ' loss = ',round(rmse),
              ' population = ',round(popCount)))+
  scale_y_continuous(limits=c(0,max_val))
  #scale_x_continuous(limits=c(earliest_obs,latest_obs))

p2 <- ggplot()
yearPred <- getPrediction(optRes$par)
max_pred_val <- max(c(obsSource$Infected,obsSource$Removed,yearPred$infected.pred,yearPred$removed.pred)) * 1.1
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