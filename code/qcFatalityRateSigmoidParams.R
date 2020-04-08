args = commandArgs(trailingOnly=TRUE)

sigmoidParamsFile <- args[1]
outputFile <- args[2]

df1 <- read.csv(sigmoidParamsFile)
isValid <- (df1$sigmoidSpaceLag >0) &
  (df1$sigmoidSpaceLag < 6) &
  (df1$fatalityRate > 0.005) &
  (df1$fatalityRate < 0.5)

k <- df1$realDaysLag / df1$sigmoidSpaceLag

print(paste0("There is ",sum(isValid)," valid fatality rate paramters out of ",nrow(df1)))

df2 <- df1[isValid,]

validSigmoidSpaceLagMed <- median(df2$sigmoidSpaceLag)
print(paste0("validSigmoidSpaceLagMed ",validSigmoidSpaceLagMed))

validFatalityRateMed <- median(df2$fatalityRate)
print(paste0("validFatalityRateMed ",validFatalityRateMed))

df1[!isValid,]$sigmoidSpaceLag <- validSigmoidSpaceLagMed
df1[!isValid,]$fatalityRate <- validFatalityRateMed
df1[!isValid,]$realDaysLag <- df1[!isValid,]$sigmoidSpaceLag * k[!isValid]

write.csv(df1,outputFile,row.names = F)
print('Done')
