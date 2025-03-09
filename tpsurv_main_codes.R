library(MASS)
library(glmnet)
library(survival)
library(DescTools)
library(mice)
library(smcfcs)
library(survcomp)
library(Hmisc)
library(knitr)

# case 1 a binary missing covariate
# case 2 a continuous missing covariate
# case 3 two missing covariates

# setting 1 MCAR
# setting 2 MAR
# setting 3 MARviol

nIter = 100
param = list(case = 1, 
             setting = 1, 
             N = 150, 
             p = 10, 
             alpha = 1.25, 
             beta = c(rep(0.5, p/2), rep(0, p/2)),
             r = 0.3,
             lambda = 1,
             rho = 1,
             c0 = 10)

res = evalRun(nIter = nIter, param = param)

modelnames = c("CCA", "NI", "MI-Wood", "MI-Bartlett", "EG")

tab = vector(mode = "list")
tab[[1]] = res$cindex
tab[[2]] = res$CS
tab[[3]] = res$IBS
tab[[4]] = res$MCC

res.avg = matrix(NA, length(modelnames), length(tab))
rownames(res.avg) = modelnames
colnames(res.avg) = c("c-index", "CS", "IBS", "MCC")

for(j in 1:length(tab)){
  res.avg[,j] = t(t(apply(tab[[j]], 1, function(x) 
    paste0(format(round(mean(x, na.rm=T),2),nsmall=2), " (", 
           format(round(sd(x, na.rm=T),2),nsmall=2), ")"))))
}

print(kable(res.avg))

