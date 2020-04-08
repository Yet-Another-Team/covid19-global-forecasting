require(pracma)
require(stringr)
require(ggplot2)

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


args = commandArgs(trailingOnly=TRUE)

sigmoidSpaceConfirmedTsFile <- args[1]
outputFile <- args[2]
outdirFigures <- args[3]

sigmoidSpaceConfirmedTsDf <- read.csv(sigmoidSpaceConfirmedTsFile)

sigmoidSpaceConfirmedTsDf$key <- getKey(sigmoidSpaceConfirmedTsDf$province, sigmoidSpaceConfirmedTsDf$country)

locationsDf <- data.frame(
  province=sigmoidSpaceConfirmedTsDf$province,
  country=sigmoidSpaceConfirmedTsDf$country)
locationsDf <- locationsDf[!duplicated(locationsDf),]

locationsDf$key <- getKey(locationsDf$province, locationsDf$country)

N <- nrow(locationsDf)


resDf <- NULL

for(i in (1:N)) {
  rowDf <- locationsDf[i,]
  province <- rowDf$province
  country <- rowDf$country
  key <- rowDf$key
  locSpecificDf <- sigmoidSpaceConfirmedTsDf[sigmoidSpaceConfirmedTsDf$key == key,]
  
  rmse <- function(obs,pred) {
    sqrt(mean((obs-pred)*(obs-pred)))
  }
  
  
  toMinimize <- function(p) {
    shift <- p[["sigmoidSpaceLag"]]
    scale <- sigmoid(p[["fatalityRate"]])
    predY <- sigmoid(locSpecificDf$x-shift)*scale
    rmse(predY, locSpecificDf$y)
  }
  
  optRes <- NULL
  
  for(j in 1:100) {
    set.seed(125433 + 101*i+31*j)
    startP = c(sigmoidSpaceLag = runif(1),fatalityRate=runif(1,min=-5,max=5)) 
    optCtr <- list(trace=0,maxit=10000)
    curOptRes <- optim(startP, toMinimize,control = optCtr)
    if(curOptRes$convergence != 0)
      next; # we analyze only converged results
    if(is.null(optRes) || optRes$value > curOptRes$value) {
      #print(paste0("Iteration ",j,": loss improved from ",optRes$value," to ",curOptRes$value))
      optRes <- curOptRes
    } else {
      #print(paste0("Iteration ",j,": loss did not improve (",curOptRes$value,")"))
    }
  }
  rmse <- optRes$value
  sigmoidSpaceLag <- optRes$par[["sigmoidSpaceLag"]]
  fatalityRate <- sigmoid(optRes$par[["fatalityRate"]])
  curResRowDf <- data.frame(province=province, country=country, loss = rmse, sigmoidSpaceLag = sigmoidSpaceLag, fatalityRate=fatalityRate, key = key)
  print(paste0("province ",province," country ",country," loss ",rmse, " sigmoidSpaceLag ",sigmoidSpaceLag,' fatalityRate ',fatalityRate))
  if(is.null(resDf)) {
    resDf <- curResRowDf
  } else {
    resDf <- rbind(resDf,curResRowDf)
  }
  
  x2 <- seq(from=-10,to=20, by=0.1)
  y2 <- sigmoid(x2 - sigmoidSpaceLag)*fatalityRate
  df2 <- data.frame(x=x2,y=y2)
  
  p <- ggplot(locSpecificDf) +
    labs(subtitle = paste0("sigmoidSpaceLag = ",round(sigmoidSpaceLag,3),
                      " fatalityRate = ",round(fatalityRate,4)))+
    geom_point(aes(x=x,y=y),shape=1) +
    geom_line(data=df2,aes(x=x,y=y))
  
  
  outFigureFile <- file.path(outdirFigures,paste0(key,".png"))
  ggsave(outFigureFile,p)
}

write.csv(resDf, outputFile,row.names = F)
print("Done")



