library(MASS)
library(glmnet)
library(survival)
library(DescTools)
library(mice)
library(smcfcs)
library(survcomp)
library(Hmisc)

# case 1 a binary missing covariate
# case 2 a continuous missing covariate
# case 3 two missing covariates

# setting 1 MCAR
# setting 2 MAR
# setting 3 MARviol

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

evalrun = function(nIter, param){
  
  nModel = 5
  cmat = csmat = ibsmat = TP = TN = FP = FN = MCC = matrix(NA, nModel, nIter)
  
  for(jj in 1:nIter){
    
    set.seed(jj)
    train = simulDat(case=param$case, setting=param$setting, N=param$N, p=param$p, alpha=param$alpha, 
                     beta=param$beta, r=param$r, lambda=param$lambda, rho=param$rho, c0=param$c0)
    # while(sum(train$targdat$status) < 7){ # to avoid an extreme case of accidentally too many censored observations left in target samples
    #   rdnum = sample(1:1000,1)
    #   set.seed(sample(jj*rdnum))
    #   train = simulDat(case=param$case, setting=param$setting, N=param$N, p=param$p, alpha=param$alpha,
    #                    beta=param$beta, r=param$r, lambda=param$lambda, rho=param$rho, c0=param$c0)
    # }
    test  = simulDat(case=param$case, setting=param$setting, N=param$N, p=param$p, alpha=param$alpha, 
                     beta=param$beta, r=param$r, lambda=param$lambda, rho=param$rho, c0=param$c0)
    
    # step 1A for expert-guided
    y = cbind(time = train$fulldat$time, status = train$fulldat$status)
    enet.ini    = cv.glmnet(train$fulldat$U, y, family = "cox", nfolds = 5, type.measure="deviance", alpha=0.5)
    abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(train$fulldat$U))
    w = abs.enet
    w[1:2] = 0
    premodel = cv.glmnet(train$fulldat$U, y, family = "cox", nfolds = 5, type.measure="deviance", penalty.factor = w)
    
    # step 1B for adding zeta
    train$targdat$zeta = predict(premodel, newx=train$targdat$U, s="lambda.1se")
    test$targdat$zeta  = predict(premodel, newx=test$targdat$U, s="lambda.1se")
    
    # setup for train data
    y = cbind(time = train$targdat$time, status = train$targdat$status)
    x = data.matrix(cbind(train$targdat$U, train$targdat$V))
    x.eg = data.matrix(cbind(train$targdat$zeta, train$targdat$V))
    
    # expert-guided
    m.eg = coxph(Surv(y) ~ x.eg)
    coef.m.eg = coef(m.eg)
    if(sum(coef(premodel, s="lambda.1se")==0)==param$p){ coef.m.eg[1] = 0}
    
    # complete-case analysis
    enet.ini    = cv.glmnet(x, y, family = "cox", nfolds = 5, type.measure="deviance", alpha=0.5)
    abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(x))
    w = abs.enet
    m.cca = cv.glmnet(x, y, family = "cox", nfolds = 5, type.measure="deviance", penalty.factor = w)
    
    # native imputation (mean)
    # new setup for train data
    y = cbind(time = train$fulldat$time, status = train$fulldat$status)
    if(param$case==1){
      imputedV.ni = as.numeric(with(train$fulldat, Hmisc::impute(V, function(x) Mode(x, na.rm=T)[1])))
    }else if(param$case==2){
      imputedV.ni = as.numeric(with(train$fulldat, Hmisc::impute(V, function(x) mean(x, na.rm=T)[1])))
    }else if(param$case==3){
      imputedV.ni1 = as.numeric(with(train$fulldat, Hmisc::impute(V[,1], function(x) Mode(x, na.rm=T)[1])))
      imputedV.ni2 = as.numeric(with(train$fulldat, Hmisc::impute(V[,2], function(x) mean(x, na.rm=T)[1])))
      imputedV.ni = cbind(imputedV.ni1, imputedV.ni2)
    }
    imputed.x.ni = data.matrix(cbind(train$fulldat$U, imputedV.ni))
    enet.ini    = cv.glmnet(imputed.x.ni, y, family = "cox", nfolds = 5, type.measure="deviance", alpha=0.5)
    abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(imputed.x.ni))
    w = abs.enet
    m.ni = cv.glmnet(imputed.x.ni, y, family = "cox", nfolds = 5, type.measure="deviance", penalty.factor = w)
    covdata = data.frame(U=train$fulldat$U, V=train$fulldat$V)
    
    # setup for test data
    testy = cbind(time = test$targdat$time, status = test$targdat$status)
    testy.mi = cbind(time = test$fulldat$time, status = test$fulldat$status)
    testx = data.matrix(cbind(test$targdat$U, test$targdat$V))
    testx.eg = data.matrix(cbind(test$targdat$zeta, test$targdat$V))
    if(param$case==1){
      imputedV.ni = as.numeric(with(test$fulldat, Hmisc::impute(V, function(x) Mode(x, na.rm=T)[1])))
    }else if(param$case==2){
      imputedV.ni = as.numeric(with(test$fulldat, Hmisc::impute(V, function(x) mean(x, na.rm=T)[1])))
    }else if(param$case==3){
      imputedV.ni1 = as.numeric(with(test$fulldat, Hmisc::impute(V[,1], function(x) Mode(x, na.rm=T)[1])))
      imputedV.ni2 = as.numeric(with(test$fulldat, Hmisc::impute(V[,2], function(x) mean(x, na.rm=T)[1])))
      imputedV.ni = cbind(imputedV.ni1, imputedV.ni2)
    }
    imputed.testx.ni = data.matrix(cbind(test$fulldat$U, imputedV.ni))
    
    # multiple imputation
    # Step 1: Impute data using mice
    M = 5
    imp = mice(covdata, m = M, printFlag = F) # uses defaultMethod (logreg for binary; pmm for continuous)
    
    # Step 2: Fit Lasso models on each imputed dataset
    mi_models = lapply(1:M, function(i) {
      x_imp = data.matrix(complete(imp, action = i))
      enet.ini    = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", alpha=0.5)
      abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(imputed.x.ni))
      w = abs.enet
      finmodel = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", penalty.factor = w) 
      return(coef(finmodel, s="lambda.1se"))
    })
    
    # Step 3: Wood et al.
    selected_wood_half = which(rowSums(do.call(cbind, mi_models)!=0)>M/2)
    mi_wood_half = lapply(1:M, function(i) {
      x_imp = data.matrix(complete(imp, action = i)[,selected_wood_half])
      if(length(selected_wood_half)>1){
        enet.ini    = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", alpha=0.5)
        abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(x_imp))
        w = abs.enet
        finmodel = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", penalty.factor = w)
        return(coef(finmodel, s="lambda.1se"))
      }else if(length(selected_wood_half)==1){
        return(coef(coxph(Surv(y) ~ x_imp)))
      }else if(length(selected_wood_half)==0){
        return(rep(0, length(selected_wood_half)))
      }
    })
    coef_mi_wood_half = rowMeans(do.call(cbind, mi_wood_half))
    
    # Step 1: Impute data using smcfcs by Bartlett et al. (2015)
    datforMI = data.frame(time = train$fulldat$time, status = train$fulldat$status, U = train$fulldat$U, V = train$fulldat$V)
    imputeindex = rep("", ncol(datforMI))
    if(param$case==1){
      imputeindex[which(colSums(is.na(datforMI))>0)] = "brlogreg"
    }else if(param$case==2){
      imputeindex[which(colSums(is.na(datforMI))>0)] = "norm"
    }else if(param$case==3){
      imputeindex[which(colSums(is.na(datforMI))>0)][1] = "brlogreg"
      imputeindex[which(colSums(is.na(datforMI))>0)][2] = "norm"
    }
    tryCatch({
      if(param$case %in% c(1:2)){
        uptMI = quiet(smcfcs(datforMI, smtype="coxph", smformula="Surv(time,status) ~ U.1+U.2+U.3+U.4+U.5+U.6+U.7+U.8+U.9+U.10+V",
                             method=imputeindex, m = 5)) 
      }else if(param$case==3){
        uptMI = quiet(smcfcs(datforMI, smtype="coxph", smformula="Surv(time,status) ~ U.1+U.2+U.3+U.4+U.5+U.6+U.7+U.8+U.9+U.10+V.1+V.2",
                             method=imputeindex, m = 5))
      }
      
      # Step 2: Fit Lasso models on each imputed dataset
      mi_models = lapply(1:M, function(i) {
        x_imp = data.matrix(uptMI$impDatasets[[i]][,-c(1:2)])
        enet.ini    = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", alpha=0.5)
        abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(imputed.x.ni))
        w = abs.enet
        finmodel = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", penalty.factor = w) 
        return(coef(finmodel, s="lambda.1se"))
      })
      
      # Step 3: Bartlett et al.
      selected_bart_half = which(rowSums(do.call(cbind, mi_models)!=0)>M/2)
      mi_bart_half = lapply(1:M, function(i) {
        x_imp = data.matrix(complete(imp, action = i)[,selected_bart_half])
        if(length(selected_bart_half)>1){
          enet.ini    = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", alpha=0.5)
          abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(x_imp))
          w = abs.enet
          finmodel = cv.glmnet(x_imp, y, family = "cox", nfolds = 5, type.measure="deviance", penalty.factor = w)
          return(coef(finmodel, s="lambda.1se"))
        }else if(length(selected_bart_half)==1){
          return(coef(coxph(Surv(y) ~ x_imp)))
        }else if(length(selected_bart_half)==0){
          return(rep(0, length(selected_bart_half)))
        }
      })
      coef_mi_bart_half = rowMeans(do.call(cbind, mi_bart_half))
    }, error=function(e){})
    
    # variable selection performance
    truecoef = c(param$beta, param$alpha)
    coefmat = matrix(0, length(truecoef), nModel)
    coefmat[,1] = as.numeric(coef(m.cca, s="lambda.1se"))
    coefmat[,2] = as.numeric(coef(m.ni, s="lambda.1se"))
    if(length(selected_wood_half)>0){coefmat[selected_wood_half,3] = coef_mi_wood_half}
    if(length(selected_bart_half)>0){coefmat[selected_bart_half,4] = coef_mi_bart_half}
    coefmat[,5] = c(as.vector(coef.m.eg[1]*coef(premodel, s="lambda.1se")), coef.m.eg[-1])
    
    TP[,jj] = colSums((truecoef!=0) & (coefmat!=0))
    TN[,jj] = colSums((truecoef==0) & (coefmat==0))
    FP[,jj] = colSums((truecoef==0) & (coefmat!=0))
    FN[,jj] = colSums((truecoef!=0) & (coefmat==0))
    MCC[,jj] = (TP[,jj]*TN[,jj] - FP[,jj]*FN[,jj]) / sqrt((TP[,jj]+FP[,jj])*(TP[,jj]+FN[,jj])*(TN[,jj]+FP[,jj])*(TN[,jj]+FN[,jj]))
    MCC[colSums(coefmat)==0,jj] = 0
    
    # compute predictions
    pred.train = pred.test = vector("list", nModel)
    
    pred.train[[1]]  = x %*% coefmat[,1] #predict(m.cca, newx = x, s="lambda.1se")
    pred.train[[2]]  = x %*% coefmat[,2]
    # if(param$case %in% c(1:2)){
    #   pred.train[[2]]  = predict(m.ni, newx = imputed.x.ni[which(!is.na(train$fulldat$V)),], s="lambda.1se")
    # }else if(param$case==3){
    #   pred.train[[2]]  = predict(m.ni, newx = imputed.x.ni[which(!is.na(train$fulldat$V[,1])),], s="lambda.1se")
    # }
    pred.train[[3]]  = x %*% coefmat[,3]
    pred.train[[4]]  = x %*% coefmat[,4]
    pred.train[[5]]  = x %*% coefmat[,5]
    
    pred.test[[1]]  = testx %*% coefmat[,1] #predict(m.cca, newx = testx, s="lambda.1se")
    pred.test[[2]]  = testx %*% coefmat[,2] #predict(m.ni, newx = testx, s="lambda.1se")
    pred.test[[3]]  = testx %*% coefmat[,3]
    pred.test[[4]]  = testx %*% coefmat[,4]
    pred.test[[5]]  = testx %*% coefmat[,5]
    
    for(k in 1:nModel){
      
      # Harrell's c-index
      cmat[k,jj] = glmnet::Cindex(pred=pred.test[[k]], y=testy)
      # Calibration slope
      csmat[k,jj] = coef(coxph( Surv(time, status) ~ pred.test[[k]], data = test$targdat))
      # Integrated brier score (ibs)
      tryCatch({
        ibsmat[k,jj] = sbrier.score2proba(data.frame(time=train$targdat$time, event=train$targdat$status, score=as.vector(pred.train[[k]])), 
                                          data.frame(time=test$targdat$time,  event=test$targdat$status,  score=as.vector(pred.test[[k]])), method="cox")$bsc.integrated
      }, error=function(e){})
    }
    print(jj)
  }
  
  return(list(cmat=cmat, csmat=csmat, ibsmat=ibsmat,
              TP=TP, TN=TN, FP=FP, FN=FN, MCC=MCC))
}