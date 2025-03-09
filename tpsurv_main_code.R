# the file is licensed under the GNU General Public License v3.0 (GPL-3.0)

# you need to set your own working directory
setwd("Specify-Your-Working-Directory")

# load packages
library(MASS)
library(glmnet)
library(survival)
library(DescTools)
library(mice)
library(smcfcs)
library(survcomp)
library(Hmisc)
library(knitr)

# parameters
# - nIter: # of Monte Carlo simulations
# - case: 1 (a binary missing covariate); 2 (a continuous missing covariate); 3 (two missing covariates)
# - setting: 1 (MCAR); 2 (MAR); 3 (MARviol)
# - n: # of total observations in two-phase data
# - p: dimension of U
# - alpha: coefficient(s) of V
# - beta: coefficients of U
# - r: ratio of two-phase data (n'/n)
# - lambda: shape parameter
# - rho: scale parameter
# - c0: parameter for censoring times 

# import several R functions (note: save these files from the folder 'functions' to your working directory)
source("simul_dat_fun")
source("run_methods_eval_fun")

# run the 
res = run_methods_eval_fun(nIter = nIter, case = case, setting = setting, n = n, p = p, 
                           alpha = alpha, beta = beta, r = r, lambda = lambda, rho = rho, c0 = c0)

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

