args = commandArgs(trailingOnly=TRUE)

firstModelParamsCSVfile <- args[1]
secondModelParamsCSVfile <- args[2]

firstDf <- read.csv(firstModelParamsCSVfile)
secondDf <- read.csv(secondModelParamsCSVfile)

m1 <- merge(firstDf,secondDf, by =c('Province','Country'), suffixes = c(".fst",".snd"))

origTTest <- t.test(m1$Loss.fst,m1$Loss.snd, paired = T)
print(origTTest)

# quality control
# we consider only runs that successfully fitted
m2 <- m1[
  m1$Gamma.fst < 0.2 &
    m1$Gamma.fst > 0.001 &
    m1$Gamma.snd < 0.2 &
    m1$Gamma.snd > 0.001
  ,]

print(paste0(nrow(m1)-nrow(m2)," rows were removed by QC. ",nrow(m2)," left"))

qcTTest <- t.test(m2$Loss.fst,m2$Loss.snd, paired = T)
print(qcTTest)