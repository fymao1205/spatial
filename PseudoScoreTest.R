
# ------------------------------------------
# Pseudo Score Test
# stage I working independence 
# ------------------------------------------

library(Rcpp)
library(Rsolnp)
library(MASS)#
#library(copula)
library(NlcOptim)
#library(nleqslv)
library(eha) #
library(parallel)
library(mnormt)#
library("pbivnorm")
sourceCpp("./sim1/commonf.cpp")
sourceCpp("./sim1/loglik.cpp")


# Pseudo Score Test 

pseudoScoreTest.f <- function(gradstep, phi1, indata_s1, XGmat, adjust.cov.Z, adjust.cov.T,
                              brks.list, marg.type="type-specific x covar-common"){
  
  len.beta <- ifelse(is.null(adjust.cov.T), 0, length(adjust.cov.T))
  
  if(marg.type=="type.specific"){
    chi.df <- len.beta*2
  }else{
    chi.df <- len.beta
  }
  
  s2 <- score_expand_phi1.f(gradstep=1e-06, phi1=c(phi1, rep(0, chi.df)), indata_s1, XGmat, 
                            adjust.cov.Z, adjust.cov.T, brks.list, marg.type) 
  nsample = length(unique(indata_s1$id))
  J2 <- (t(s2) %*% s2)/nsample
  
  est.pst <- (colSums(s2) %*% ginv(J2) %*% colSums(s2))/nsample
  
  pval=1-pchisq(est.pst, df=chi.df)
  
  res.list <- list(pst.stat=est.pst,
                   pval=pval)
  return(res.list)
}

eval.obs.logL_expand_j.f <- function(l, r, cov, alp, pi_marg, brks){
  
  comp.L_1j <- Lcomp_adjusted_1j(l, r, cov, alp, brks)
  use.finite <- is.finite(r)
  pi_marg_0 <- 1- pi_marg
  if(any(use.finite)){ pi_marg_0[use.finite] <- 0 }
  
  res <- log(comp.L_1j*pi_marg + pi_marg_0)
  
  res <- ifelse(is.nan(res), .Machine$double.xmin, res)
  return(res)
} 

lik_expand_phi1.f <- function(phi1, indata_s1, XGmat, adjust.cov.Z, adjust.cov.T, brks.list, 
                              marg.type="type-specific x covar-common"){
  
  # J=2; For J>2, expand in a similar spirit...
  indata_s1_j1 = subset(indata_s1, type==1)
  indata_s1_j2 = subset(indata_s1, type==2)
  
  R1 = length(brks.list[[1]]) +1; R2 = length(brks.list[[2]]) +1; 
  R = R1 + R2 
  
  XGmat_j1 <- as.matrix(subset(XGmat, indata_s1$type==1))
  XGmat_j2 <- as.matrix(subset(XGmat, indata_s1$type==2))
  
  len.eta1 <- (!is.null(adjust.cov.Z))*length(adjust.cov.Z)
  len.beta <- (!is.null(adjust.cov.T))*length(adjust.cov.T)
  
  if(marg.type=="common"){
    
    alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[(1:R1)]); 
    eta_10 = phi1[R1+1]; eta_20 = phi1[R1+1];  
    
    if(len.eta1>0){
      eta_21 = eta_11 = phi1[c(R1+1+(1:len.eta1))]
      pi_marg_1 = expit.f(eta_10+XGmat_j1[,adjust.cov.Z, drop=F] %*% eta_11)
      pi_marg_2 = expit.f(eta_20+XGmat_j2[,adjust.cov.Z, drop=F] %*% eta_21)
    }else{
      pi_marg_1 = expit.f(rep(eta_10, dim(XGmat_j1)[1]))
      pi_marg_2 = expit.f(rep(eta_20, dim(XGmat_j2)[1]))
    }
    
    if(len.beta>0){
      beta_21 = beta_11 = phi1[c(R1+1+(len.eta1)+(1:len.beta))]
      cov_1 = XGmat_j1[,adjust.cov.T,drop=F] %*% beta_11
      cov_2 = XGmat_j2[,adjust.cov.T,drop=F] %*% beta_21
    }else{
      cov_1 = rep(0.0, dim(XGmat_j1)[1])
      cov_2 = rep(0.0, dim(XGmat_j2)[1])
    }
    
    
    obs.logL_1 <- eval.obs.logL_expand_j.f(indata_s1_j1$l, indata_s1_j1$r, cov_1, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_expand_j.f(indata_s1_j2$l, indata_s1_j2$r, cov_2, alp_2, pi_marg_2, brks.list[[2]])
    
  }else if(marg.type=="type-specific x covar-common"){
    
    alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[R1+1+(1:R2)]); 
    eta_10 = phi1[R1+1]; eta_20 = phi1[R1+R2+2];  
    
    if(len.eta1>0){
      eta_21 = eta_11 = phi1[c(R1+R2+2+(1:len.eta1))]
      pi_marg_1 = expit.f(eta_10+XGmat_j1[,adjust.cov.Z, drop=F] %*% eta_11)
      pi_marg_2 = expit.f(eta_20+XGmat_j2[,adjust.cov.Z, drop=F] %*% eta_21)
    }else{
      pi_marg_1 = expit.f(rep(eta_10, dim(XGmat_j1)[1]))
      pi_marg_2 = expit.f(rep(eta_20, dim(XGmat_j2)[1]))
    }
    
    if(len.beta>0){
      beta_21 = beta_11 = phi1[c(R1+R2+2+(len.eta1)+(1:len.beta))]
      cov_1 = XGmat_j1[,adjust.cov.T,drop=F] %*% beta_11
      cov_2 = XGmat_j2[,adjust.cov.T,drop=F] %*% beta_21
    }else{
      cov_1 = rep(0.0, dim(XGmat_j1)[1])
      cov_2 = rep(0.0, dim(XGmat_j2)[1])
    }
    
    obs.logL_1 <- eval.obs.logL_expand_j.f(indata_s1_j1$l, indata_s1_j1$r, cov_1, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_expand_j.f(indata_s1_j2$l, indata_s1_j2$r, cov_2, alp_2, pi_marg_2, brks.list[[2]])
    
  }else if(marg.type=="type-specific"){
    
    alp_1 = exp(phi1[1:R1]);  
    eta_10 = phi1[R1+1];   
    #eta_11 = phi1[c(R1+1+(1:ncol_Gmat))]; 
    alp_2 = exp(phi1[R1+1+len.eta1+(1:R2)]); 
    eta_20 = phi1[R1+R2+2+len.eta1];
    
    if(len.eta1>0){
      eta_11 = phi1[c(R1+1+(1:len.eta1))]
      eta_21 = eta_11 = phi1[c(R1+R2+2+len.eta1+(1:len.eta1))]
      pi_marg_1 = expit.f(eta_10+XGmat_j1[,adjust.cov.Z, drop=F] %*% eta_11)
      pi_marg_2 = expit.f(eta_20+XGmat_j2[,adjust.cov.Z, drop=F] %*% eta_21)
    }else{
      pi_marg_1 = expit.f(rep(eta_10, dim(XGmat_j1)[1]))
      pi_marg_2 = expit.f(rep(eta_20, dim(XGmat_j2)[1]))
    }
    
    if(len.beta>0){
      beta_11 = phi1[c(R1+R2+2+(len.eta1)*2+(1:len.beta))]
      beta_21 = phi1[c(R1+R2+2+(len.eta1)*2+len.beta+(1:len.beta))]
      cov_1 = XGmat_j1[,adjust.cov.T,drop=F] %*% beta_11
      cov_2 = XGmat_j2[,adjust.cov.T,drop=F] %*% beta_21
    }else{
      cov_1 = rep(0.0, dim(XGmat_j1)[1])
      cov_2 = rep(0.0, dim(XGmat_j2)[1])
    }
    
    obs.logL_1 <- eval.obs.logL_expand_j.f(indata_s1_j1$l, indata_s1_j1$r, cov_1, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_expand_j.f(indata_s1_j2$l, indata_s1_j2$r, cov_2, alp_2, pi_marg_2, brks.list[[2]])
  }
  
  res <- c(obs.logL_1, obs.logL_2)
  
  re_indata_s1 <- as.data.frame(rbind(indata_s1_j1, indata_s1_j2))
  
  res <- aggregate(res, list(re_indata_s1$id), sum, drop = T, simplify = T)[, -1]
  
  return(res)
  
}

score_expand_phi1.f <- function(gradstep, phi1, indata_s1, XGmat, adjust.cov.Z, adjust.cov.T,
                                brks.list, marg.type="type-specific x covar-common"){
  
  lik_phi1 <- lik_expand_phi1.f(phi1, indata_s1, XGmat, adjust.cov.Z, adjust.cov.T, brks.list, marg.type)
  
  p1 = length(phi1)
  #score_mat <- matrix(NA, nrow=nsample, ncol=p1)
  
  score.list <- lapply(1:p1, function(x){
    
    phi1p = phi1
    phi1p[x] = phi1[x] + gradstep
    
    profile_x <- lik_expand_phi1.f(phi1p, indata_s1, XGmat, adjust.cov.Z, adjust.cov.T, brks.list, marg.type)
    tmp <- (profile_x - lik_phi1)/gradstep
    
    tmp
  })
  
  score_mat <- do.call("cbind", score.list)
  
  return(score_mat)
}

est_s1_H0.f <- function(indata_s1, XGmat, adjust.cov.Z, brks.list, marg.type="type-specific x covar-common"){
  
  if(is.null(marg.type)){
    marg.type = "type-specific x covar-common"
  }

  len.eta1 <- (!is.null(adjust.cov.Z))*length(adjust.cov.Z)
  
  J = max(indata_s1$type)
  indata_s1.list = lapply(1:J, function(j){
    subset(indata_s1, type==j)
  })
  
  XGmat.list = lapply(1:J, function(j){
    as.matrix(subset(XGmat, indata_s1$type==j))
  })
  
  R_vec <- sapply(brks.list, length)+1

  if(marg.type == "type-specific x covar-common"){
    
    obj.f <- function(phi1){
      
      eta_1 = phi1[sum(R_vec+1)+1:len.eta1]
      
      res1 <- lapply(1:J, function(j){
        
        alp =  exp(phi1[c(0,cumsum(R_vec+1))[j]+(1:R_vec[j])]); 
        eta_0 = phi1[cumsum(R_vec+1)[j]]
        xgmat <- XGmat.list[[j]]
        pi_marg = expit.f(eta_0+ xgmat[,adjust.cov.Z, drop=F] %*% eta_1)
        dt=indata_s1.list[[j]]
        obs.logL <- eval.obs.logL_j.f(dt$l, dt$r, alp, pi_marg, brks.list[[j]])
        res0 = -sum(obs.logL)
        #print(res0)
        res0
      })
      
      res <- do.call("sum", res1)
      
      res <- ifelse(is.finite(res), res, .Machine$double.xmax)
      
      
      print(phi1); print(res)
      
      return(res)
    }
    
    dim_phi1 = sum(R_vec+1)+len.eta1
    system.time(
      #res.s1 <- nlm(f=obj.f, p=rep(-1, dim_phi1),gradtol=1e-04, steptol=1e-05, hessian=TRUE)
      res.s1 <- optim(par=rep(-1, dim_phi1), fn=obj.f, 
                      control = list(maxit= 1000), method="L-BFGS-B", hessian = TRUE)
    )
    
  }
  
  if(marg.type == "common"){
    
    R <- length(brks.list[[1]])+1
   
    obj.f <- function(phi1){
      
      alp = exp(phi1[1:R]);
      eta = phi1[-(1:R)]; 
      #pi_marg = expit.f(eta[1]+eta[2]*indata_s1$x+eta[3]*indata_s1$g)
      pi_marg = expit.f(eta[1]+ XGmat[,adjust.cov.Z, drop=F] %*% eta[-1])
      
      obs.logL <- eval.obs.logL_j.f(indata_s1$l, indata_s1$r, alp, pi_marg, brks.list[[1]])
      
      res <- -sum(obs.logL)
      
      print(res)
      
      return(res)
    }
    
    dim_phi1 <- R+1+len.eta1
    res.s1 <- optim(par=rep(-1, dim_phi1), fn=obj.f,control = list(maxit= 1000),  method="L-BFGS-B", hessian=TRUE)
  }
  
  if(marg.type == "type-specific"){
    
    obj.f <- function(phi1){
      
      res1 <- lapply(1:J, function(j){
        
        #alp =  exp(phi1[c(0,R_vec)[j]+(1:R_vec[j])]); 
        exp(phi1[c(0,cumsum(R_vec+1))[j]+(1:R_vec[j])]); 
        eta_0 = phi1[cumsum(R_vec+1)[j]]
        eta_1 = phi1[sum(R_vec+1)+(j-1)*len.eta1+1:len.eta1]
        xgmat <- XGmat.list[[j]]
        pi_marg = expit.f(eta_0+ xgmat[,adjust.cov.Z, drop=F] %*% eta_1)
        dt=indata_s1.list[[j]]
        obs.logL <- eval.obs.logL_j.f(dt$l, dt$r, alp, pi_marg, brks.list[[j]])
        res0 = -sum(obs.logL)
        res0
      })
      
      res <- do.call("sum", res1)
      
      res <- ifelse(is.finite(res), res, .Machine$double.xmax)
      return(res)
      
      return(res)
    }
    
    dim_phi1 = sum(R_vec+1+len.eta1)
    system.time(
      #res.s1 <- nlm(f=obj.f, p=rep(-1, dim_phi1),gradtol=1e-04, steptol=1e-05, hessian=TRUE)
      res.s1 <- optim(par=rep(-1, dim_phi1), fn=obj.f, 
                      control = list(maxit= 1000), method="L-BFGS-B", hessian = TRUE)
    )
    
  }
  
  return(res.s1)
}



# test
#indata_s1 = dt_s1
#adjust.cov.Z = 1:length(covar_str)
#XGmat <- dt_s1[,c(paste0("x", adjust.cov.Z))]
#marg.type="type-specific x covar-common"
#phi.est <- est_s1_H0.f(indata_s1, XGmat, adjust.cov.Z, brks.list, marg.type="type-specific x covar-common")
#adjust.cov.T <- c(1,2, 3:length(covar_str))
#pseudoScoreTest.f(1e-06, phi.est$par, indata_s1, XGmat, adjust.cov.Z, adjust.cov.T, brks.list, marg.type="type-specific x covar-common")

#len.beta <- ifelse(is.null(adjust.cov.T), 0, length(adjust.cov.T))
#s2 <- score_expand_phi1.f(gradstep=1e-06, phi1=c(phi.est$par, rep(0, len.beta)), indata_s1, XGmat, 
#                          adjust.cov.Z, adjust.cov.T, brks.list, marg.type)
#nsample = length(unique(indata_s1$id))
#J2 <- (t(s2) %*% s2)/nsample
#est.pst <- (colSums(s2) %*% ginv(J2) %*% colSums(s2))/nsample
#pval=1-pchisq(est.pst, df=len.beta)
#pval










