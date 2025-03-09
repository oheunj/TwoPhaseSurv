simul_dat_fun = function(case, setting, n, p, alpha, beta, r, lambda, rho, c0){
  
  # p-dimensional standard normal random vectors
  U = mvrnorm(n, mu=rep(0, times=p), Sigma = diag(rep(1, p)))
  
  if(case==1){
    
    # a binary missing covariate
    prob_V = rep(NA, n)
    prob_V[which(U[,1]< -0.5)] = 0.3
    prob_V[which(U[,1]>= -0.5 & U[,1]< 1)] = 0.5
    prob_V[which(U[,1]>= 1)] = 0.2
    V = ifelse(runif(n) <= prob_V, 1, 0)
    
  }else if(case==2){
    
    # a continuous missing covariate
    V = rnorm(n, 0.4*abs(U[,1])-0.1, 0.2^2)
    
  }else if(case==3){
    
    # two missing covariates
    prob_V = rep(NA, n)
    prob_V[which(U[,1]< -0.5)] = 0.3
    prob_V[which(U[,1]>= -0.5 & U[,1]< 1)] = 0.5
    prob_V[which(U[,1]>= 1)] = 0.2
    V1 = ifelse(runif(n) <= prob_V, 1, 0)
    V2 = rnorm(0.4*abs(U[,1])-0.1, 0.2^2)
    V = cbind(V1, V2)
    colnames(V) = NULL
    
  }
  
  # Weibull latent event times
  lambda_wiki = lambda^(-1/rho)
  if(case %in% c(1,2)){
    lambda_prime = lambda^(-1/rho) / exp((U%*%beta + V*alpha)/rho)
  }else if(case==3){
    lambda_prime = lambda^(-1/rho) / exp((U%*%beta + V%*%alpha)/rho)
  }
  
  # follow-up times and event indicators
  Tlat = rweibull(n, shape=rho, scale=lambda_prime)
  C = rexp(n=n, rate=c0)
  time = pmin(Tlat, C)
  status = as.numeric(Tlat <= C)
  
  if(setting==1){
    
    prob_missing = 1-r
    missing_mask = which(runif(n) < prob_missing)
    
  }else if(setting==2){
    
    prob_missing = rep(NA, n)
    prob_missing[which(U[,1]> qnorm(r/3))]  = (1-r)/pnorm(-qnorm(r/3))
    prob_missing[which(U[,1]<= qnorm(r/3))] = 0
    missing_mask = which(runif(n) < prob_missing)
    
  }else if(setting==3){
    
    prob_missing = rep(NA, n)
    prob_missing[which(U[,1]> qnorm(r/3))]  = 
      (1-as.numeric(V[which(U[,1]> qnorm(r/3))] <= 0)*0.1-r)/pnorm(-qnorm(r/3)) + 
      as.numeric(V[which(U[,1]> qnorm(r/3))] <= 0)*0.1
    prob_missing[which(U[,1]<= qnorm(r/3))] = 
      as.numeric(V[which(U[,1]<= qnorm(r/3))] <= 0)*0.1
    missing_mask = which(runif(n) < prob_missing)
    
  }
  
  targU = U[-missing_mask,]
  if(case %in% c(1,2)){
    targV = V[-missing_mask]
    V[missing_mask] = NA
  }else if(case==3){
    targV = V[-missing_mask,]
    V[missing_mask,] = rep(NA, ncol(V))
  }
  targtime = time[-missing_mask]
  targstatus = status[-missing_mask]
  
  fulldat = list(time=time, status=status, U=U, V=V)
  targdat = list(time=targtime, status=targstatus, U=targU, V=targV)
  
  return(list(fulldat=fulldat, targdat=targdat))
}
