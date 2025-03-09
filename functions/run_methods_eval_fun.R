run_methods_eval_fun = function(nIter, param){
  # to disable in-function printed message for smcfcs
  quiet = function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }
  nModel = 5
  cindex = CS = IBS = MCC = matrix(NA, nModel, nIter)
  
  for(jj in 1:nIter){
    
    set.seed(jj)
    train = simulDat(case = param$case, setting = param$setting, N = param$N, p = param$p, alpha = param$alpha, 
                     beta = param$beta, r = param$r, lambda = param$lambda, rho = param$rho, c0 = param$c0)
    while(sum(train$targdat$status) < 7){ # avoid an extreme case of accidentally too many censored observations left in target samples
      rdnum = sample(1:1000,1)
      set.seed(sample(jj*rdnum))
      train = simulDat(case=param$case, setting=param$setting, N=param$N, p=param$p, alpha=param$alpha,
                       beta=param$beta, r=param$r, lambda=param$lambda, rho=param$rho, c0=param$c0)
    }
    test  = simulDat(case = param$case, setting = param$setting, N = param$N, p = param$p, alpha = param$alpha, 
                     beta = param$beta, r = param$r, lambda = param$lambda, rho = param$rho, c0 = param$c0)
    
    # step 1A for expert-guided
    y = cbind(time = train$fulldat$time, status = train$fulldat$status)
    U = train$fulldat$U
    enet.ini = cv.glmnet(U, y, family = "cox", nfolds = 5, type.measure = "deviance", alpha = 0.5)
    abs.enet = 1/(abs(coef(enet.ini)) + 1/nrow(U))
    w = abs.enet
    w[1:2] = 0
    premodel = cv.glmnet(U, y, family = "cox", nfolds = 5, type.measure = "deviance", penalty.factor = w)
    
    # step 1B for adding zeta
    train$targdat$zeta = predict(premodel, newx = train$targdat$U, s = "lambda.1se")
    test$targdat$zeta  = predict(premodel, newx = test$targdat$U, s = "lambda.1se")
    
    # setup for target samples using train data
    y = cbind(time = train$targdat$time, status = train$targdat$status)
    X = data.matrix(cbind(train$targdat$U, train$targdat$V))
    X.eg = data.matrix(cbind(train$targdat$zeta, train$targdat$V))
    
    # expert-guided
    m.eg = coxph(Surv(y) ~ X.eg)
    coef.m.eg = coef(m.eg)
    if(sum(coef(premodel, s = "lambda.1se")==0)==param$p){ coef.m.eg[1] = 0 }
    
    # complete-case analysis
    enet.ini = cv.glmnet(X, y, family = "cox", nfolds = 5, type.measure = "deviance", alpha = 0.5)
    abs.enet = 1/(abs(coef(enet.ini)) + 1/nrow(X))
    w = abs.enet
    m.cca = cv.glmnet(X, y, family = "cox", nfolds = 5, type.measure = "deviance", penalty.factor = w)
    
    # native imputation (mean)
    y = cbind(time = train$fulldat$time, status = train$fulldat$status)
    if(param$case==1){
      imputed.V.ni = as.numeric(with(train$fulldat, Hmisc::impute(V, function(x) Mode(x, na.rm=T)[1])))
    }else if(param$case==2){
      imputed.V.ni = as.numeric(with(train$fulldat, Hmisc::impute(V, function(x) mean(x, na.rm=T)[1])))
    }else if(param$case==3){
      imputed.V.ni1 = as.numeric(with(train$fulldat, Hmisc::impute(V[,1], function(x) Mode(x, na.rm=T)[1])))
      imputed.V.ni2 = as.numeric(with(train$fulldat, Hmisc::impute(V[,2], function(x) mean(x, na.rm=T)[1])))
      imputed.V.ni  = cbind(imputed.V.ni1, imputed.V.ni2)
    }
    imputed.X.ni = data.matrix(cbind(train$fulldat$U, imputed.V.ni))
    
    enet.ini = cv.glmnet(imputed.X.ni, y, family = "cox", nfolds = 5, type.measure = "deviance", alpha = 0.5)
    abs.enet = 1/(abs(coef(enet.ini)) + 1/nrow(imputed.X.ni))
    w = abs.enet
    m.ni = cv.glmnet(imputed.X.ni, y, family = "cox", nfolds = 5, type.measure = "deviance", penalty.factor = w)
    
    # setup for test data
    testy = cbind(time = test$targdat$time, status = test$targdat$status)
    testy.mi = cbind(time = test$fulldat$time, status = test$fulldat$status)
    testX = data.matrix(cbind(test$targdat$U, test$targdat$V))
    testX.eg = data.matrix(cbind(test$targdat$zeta, test$targdat$V))
    if(param$case==1){
      imputed.V.ni = as.numeric(with(test$fulldat, Hmisc::impute(V, function(x) Mode(x, na.rm=T)[1])))
    }else if(param$case==2){
      imputed.V.ni = as.numeric(with(test$fulldat, Hmisc::impute(V, function(x) mean(x, na.rm=T)[1])))
    }else if(param$case==3){
      imputed.V.ni1 = as.numeric(with(test$fulldat, Hmisc::impute(V[,1], function(x) Mode(x, na.rm=T)[1])))
      imputed.V.ni2 = as.numeric(with(test$fulldat, Hmisc::impute(V[,2], function(x) mean(x, na.rm=T)[1])))
      imputed.V.ni = cbind(imputed.V.ni1, imputed.V.ni2)
    }
    imputed.testX.ni = data.matrix(cbind(test$fulldat$U, imputed.V.ni))
    
    # MI-Wood
    # step 1 for imputing data using mice()
    M = 5
    covdata = data.frame(U = train$fulldat$U, V = train$fulldat$V)
    imp = mice(covdata, m = M, printFlag = F)
    
    # step 2 for fitting models on each imputed dataset
    mi_models = lapply(1:M, function(i) {
      X_imp = data.matrix(complete(imp, action = i))
      enet.ini = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", alpha = 0.5)
      abs.enet = 1/(abs(coef(enet.ini)) + 1/nrow(X_imp))
      w = abs.enet
      finmodel = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", penalty.factor = w) 
      return(coef(finmodel, s = "lambda.1se"))
    })
    
    # step 3 based on a majority rule
    selected_wood = which(rowSums(do.call(cbind, mi_models)!=0)>M/2)
    mi_wood = lapply(1:M, function(i) {
      X_imp = data.matrix(complete(imp, action = i)[,selected_wood])
      if(length(selected_wood)>1){
        enet.ini = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", alpha = 0.5)
        abs.enet = 1/(abs(coef(enet.ini)) + 1/nrow(X_imp))
        w = abs.enet
        finmodel = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", penalty.factor = w)
        return(coef(finmodel, s = "lambda.1se"))
      }else if(length(selected_wood)==1){
        return(coef(coxph(Surv(y) ~ X_imp)))
      }else if(length(selected_wood)==0){
        return(rep(0, length(selected_wood)))
      }
    })
    coef_mi_wood = rowMeans(do.call(cbind, mi_wood))
    
    # MI-Bartlett
    # step 1 for imputing data using smcfcs()
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
        uptMI = quiet(smcfcs(datforMI, smtype = "coxph", 
                             smformula = paste0("Surv(time,status)~", paste(paste0("U.",1:p), collapse = "+"), "+V"),
                             method=imputeindex, m = 5)) 
      }else if(param$case==3){
        uptMI = quiet(smcfcs(datforMI, smtype = "coxph", 
                             smformula = paste0("Surv(time,status)~", paste(paste0("U.",1:p), collapse = "+"), "+V.1+V.2"),
                             method=imputeindex, m = 5))
      }
      
      # step 2 for fitting models on each imputed dataset
      mi_models = lapply(1:M, function(i) {
        X_imp = data.matrix(uptMI$impDatasets[[i]][,-c(1:2)])
        enet.ini = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", alpha = 0.5)
        abs.enet = 1/(abs(coef(enet.ini)) + 1/nrow(X_imp))
        w = abs.enet
        finmodel = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", penalty.factor = w) 
        return(coef(finmodel, s = "lambda.1se"))
      })
      
      # step 3 based on a majority rule
      selected_bart = which(rowSums(do.call(cbind, mi_models)!=0)>M/2)
      mi_bart = lapply(1:M, function(i) {
        X_imp = data.matrix(complete(imp, action = i)[,selected_bart])
        if(length(selected_bart)>1){
          enet.ini    = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", alpha = 0.5)
          abs.enet    = 1/(abs(coef(enet.ini)) + 1/nrow(X_imp))
          w = abs.enet
          finmodel = cv.glmnet(X_imp, y, family = "cox", nfolds = 5, type.measure = "deviance", penalty.factor = w)
          return(coef(finmodel, s = "lambda.1se"))
        }else if(length(selected_bart)==1){
          return(coef(coxph(Surv(y) ~ X_imp)))
        }else if(length(selected_bart)==0){
          return(rep(0, length(selected_bart)))
        }
      })
      coef_mi_bart = rowMeans(do.call(cbind, mi_bart))
    }, error=function(e){})
    
    # variable selection performance
    truecoef = c(param$beta, param$alpha)
    coefmat = matrix(0, length(truecoef), nModel)
    coefmat[,1] = as.numeric(coef(m.cca, s = "lambda.1se"))
    coefmat[,2] = as.numeric(coef(m.ni, s = "lambda.1se"))
    if(length(selected_wood)>0){coefmat[selected_wood,3] = coef_mi_wood}
    if(length(selected_bart)>0){coefmat[selected_bart,4] = coef_mi_bart}
    coefmat[,5] = c(as.vector(coef.m.eg[1]*coef(premodel, s = "lambda.1se")), coef.m.eg[-1])
    
    TP = colSums((truecoef != 0) & (coefmat != 0))
    TN = colSums((truecoef == 0) & (coefmat == 0))
    FP = colSums((truecoef == 0) & (coefmat != 0))
    FN = colSums((truecoef != 0) & (coefmat == 0))
    MCC[,jj] = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    MCC[colSums(coefmat)==0,jj] = 0
    
    # compute predictions
    pred.train = pred.test = vector("list", nModel)
    
    pred.train[[1]] = X %*% coefmat[,1]
    pred.train[[2]] = X %*% coefmat[,2]
    pred.train[[3]] = X %*% coefmat[,3]
    pred.train[[4]] = X %*% coefmat[,4]
    pred.train[[5]] = X %*% coefmat[,5]
    
    pred.test[[1]] = testX %*% coefmat[,1]
    pred.test[[2]] = testX %*% coefmat[,2]
    pred.test[[3]] = testX %*% coefmat[,3]
    pred.test[[4]] = testX %*% coefmat[,4]
    pred.test[[5]] = testX %*% coefmat[,5]
    
    for(k in 1:nModel){
      
      # Harrell's c-index
      cindex[k,jj] = glmnet::Cindex(pred=pred.test[[k]], y=testy)
      
      # Calibration slope
      CS[k,jj] = coef(coxph( Surv(time, status) ~ pred.test[[k]], data = test$targdat))
      
      # Integrated brier score (ibs)
      tryCatch({
        IBS[k,jj] = sbrier.score2proba(data.frame(time = train$targdat$time, event = train$targdat$status, 
                                                  score = as.vector(pred.train[[k]])), 
                                       data.frame(time = test$targdat$time, event = test$targdat$status,  
                                                  score = as.vector(pred.test[[k]])), 
                                       method = "cox")$bsc.integrated*10
      }, error=function(e){})
    }
    print(jj)
  }
  
  return(list(cindex = cindex, CS = CS, IBS = IBS, TP = TP, TN = TN, FP = FP, FN = FN, MCC = MCC))
}
