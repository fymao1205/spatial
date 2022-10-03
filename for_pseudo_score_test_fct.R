
eval.obs.logL_adjusted_j.f <- function(l, r, cov, alp, pi_marg, brks){
  
  comp.L_1j <- Lcomp_adjusted_1j(l, r, cov, alp, brks)
  use.finite <- is.finite(r)
  pi_marg_0 <- 1- pi_marg
  if(any(use.finite)){ pi_marg_0[use.finite] <- 0 }
  
  res <- log(comp.L_1j*pi_marg + pi_marg_0)
  
  res <- ifelse(is.nan(res), .Machine$double.xmin, res)
  return(res)
}



realdt_lik_adjusted_phi1.f <- function(phi1, indata_s1, adjust.cov.Z, adjust.cov.T, brks.list, 
                                       marg.type="type-specific x covar-common"){
  
  indata_s1_j1 = subset(indata_s1, type==1)
  indata_s1_j2 = subset(indata_s1, type==2)
  #indata_s1_j3 = subset(indata_s1, type==3)
  
  R1 = length(brks.list[[1]]) +1; R2 = length(brks.list[[2]]) +1; #R3 = length(brks.list[[3]]) +1
  R = R1 + R2 #+ R3
  
  ncol_Gmat <- dim(indata_s1_j1[,-(1:7)])[2]
  XGmat_j1 <- as.matrix(indata_s1_j1[,-(1:7)])
  XGmat_j2 <- as.matrix(indata_s1_j2[,-(1:7)])
  #XGmat_j3 <- as.matrix(indata_s1_j3[,-(1:7)])
  
  if(marg.type=="common"){
    
    alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[(1:R1)]); #alp_3 = exp(phi1[R1+R2+2+(1:R3)])
    eta_10 = phi1[R1+1]; eta_20 = phi1[R1+1];  #eta_30 = phi1[R+3]
    #eta_31 = 
    eta_21 = eta_11 = phi1[c(R1+1+(1:ncol_Gmat))]
    pi_marg_1 = expit.f(eta_10+XGmat_j1[,adjust.cov.Z, drop=F] %*% eta_11)
    pi_marg_2 = expit.f(eta_20+XGmat_j2[,adjust.cov.Z, drop=F] %*% eta_21)
    #pi_marg_3 = expit.f(eta_30+XGmat_j3 %*% eta_31)
    
    beta_21 = beta_11 = phi1[c(R1+1+(ncol_Gmat)+(1:length(adjust.cov.T)))]
    cov_1 = XGmat_j1[,adjust.cov.T,drop=F] %*% beta_11
    cov_2 = XGmat_j2[,adjust.cov.T,drop=F] %*% beta_21
    
    obs.logL_1 <- eval.obs.logL_adjusted_j.f(indata_s1_j1$l, indata_s1_j1$r, cov_1, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_adjusted_j.f(indata_s1_j2$l, indata_s1_j2$r, cov_2, alp_2, pi_marg_2, brks.list[[2]])
    #obs.logL_3 <- eval.obs.logL_j.f(indata_s1_j3$l, indata_s1_j3$r, alp_3, pi_marg_3, brks.list[[3]])
    
  }else if(marg.type=="type-specific x covar-common"){
    
    alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[R1+1+(1:R2)]); #alp_3 = exp(phi1[R1+R2+2+(1:R3)])
    eta_10 = phi1[R1+1]; eta_20 = phi1[R1+R2+2];  #eta_30 = phi1[R+3]
    #eta_31 = 
    eta_21 = eta_11 = phi1[c(R+2+(1:length(adjust.cov.Z)))]
    pi_marg_1 = expit.f(eta_10+XGmat_j1[,adjust.cov.Z, drop=F] %*% eta_11)
    pi_marg_2 = expit.f(eta_20+XGmat_j2[,adjust.cov.Z, drop=F] %*% eta_21)
    #pi_marg_3 = expit.f(eta_30+XGmat_j3 %*% eta_31)
    
    beta_21 = beta_11 = phi1[c(R+2+length(adjust.cov.Z)+(1:length(adjust.cov.T)))]
    cov_1 = XGmat_j1[,adjust.cov.T,drop=F] %*% beta_11
    cov_2 = XGmat_j2[,adjust.cov.T,drop=F] %*% beta_21
    
    obs.logL_1 <- eval.obs.logL_adjusted_j.f(indata_s1_j1$l, indata_s1_j1$r, cov_1, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_adjusted_j.f(indata_s1_j2$l, indata_s1_j2$r, cov_2, alp_2, pi_marg_2, brks.list[[2]])
    #obs.logL_3 <- eval.obs.logL_j.f(indata_s1_j3$l, indata_s1_j3$r, alp_3, pi_marg_3, brks.list[[3]])
    
  }else if(marg.type=="type-specific"){
    
    alp_1 = exp(phi1[1:R1]);  #alp_3 = exp(phi1[R1+R2+2+(1:R3)])
    eta_10 = phi1[R1+1];   #eta_30 = phi1[R+3]
    eta_11 = phi1[c(R1+1+(1:ncol_Gmat))]; 
    alp_2 = exp(phi1[R1+1+ncol_Gmat+(1:R2)]); 
    eta_20 = phi1[R1+R2+ncol_Gmat+2];
    eta_21 = phi1[c(R+2+ncol_Gmat+(1:ncol_Gmat))]
    #eta_31 = phi1[c(R+4, R+4+2*ncol_Gmat+(1:ncol_Gmat))]
    pi_marg_1 = expit.f(eta_10+XGmat_j1[,adjust.cov.Z, drop=F] %*% eta_11)
    pi_marg_2 = expit.f(eta_20+XGmat_j2[,adjust.cov.Z, drop=F] %*% eta_21)
    #pi_marg_3 = expit.f(eta_30+XGmat %*% eta_31)
    
    beta_11 = phi1[c(R+2+(ncol_Gmat)*2+(1:length(adjust.cov.T)))]
    beta_21 = phi1[c(R+2+(ncol_Gmat)*2+length(adjust.cov.T)+(1:length(adjust.cov.T)))]
    cov_1 = XGmat_j1[,adjust.cov.T,drop=F] %*% beta_11
    cov_2 = XGmat_j2[,adjust.cov.T,drop=F] %*% beta_21
    
    obs.logL_1 <- eval.obs.logL_adjusted_j.f(indata_s1_j1$l, indata_s1_j1$r, cov_1, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_adjusted_j.f(indata_s1_j2$l, indata_s1_j2$r, cov_2, alp_2, pi_marg_2, brks.list[[2]])
    #obs.logL_3 <- eval.obs.logL_j.f(indata_s1_j3$l, indata_s1_j3$r, alp_3, pi_marg_3, brks.list[[3]])
    
  }
  
  res <- c(obs.logL_1, obs.logL_2)
  
  re_indata_s1 <- as.data.frame(rbind(indata_s1_j1, indata_s1_j2))
  
  res <- aggregate(res, list(re_indata_s1$id), sum, drop = T, simplify = T)[, -1]
  
  return(res)
  
}



realdt_est_s1_adjusted.f <- function(indata_s1.list, adjust.cov.Z, adjust.cov.T,
                                     brks.list, marg.type=NULL){
  
  indata_s1 = do.call("rbind", indata_s1.list)
  J = max(indata_s1$type)
  
  R_vec <- sapply(brks.list, length)+1
  R=R_vec[1]
  
  len.eta1 <- (!is.null(adjust.cov.Z))*length(adjust.cov.Z)
  len.beta <- (!is.null(adjust.cov.T))*length(adjust.cov.T)

  if(is.null(marg.type)){
    marg.type = "common"
  }
  
  if(marg.type == "common"){
    
    R <- length(brks.list[[1]])+1
    
    obj.f <- function(phi1){
      
      alp = exp(phi1[1:R]);
      eta = phi1[R+1+(1:len.eta1)]; 
      beta = phi1[R+1+len.eta1+(1:len.beta)]
      #pi_marg = expit.f(eta[1]+eta[2]*indata_s1$x+eta[3]*indata_s1$g)
      pi_marg = expit.f(eta[1]+ XGmat[,adjust.cov.Z, drop=F] %*% eta[-1])
      cov = XGmat[,adjust.cov.T,drop=F] %*% beta
      
      obs.logL <- eval.obs.logL_adjusted_j.f(indata_s1$l, indata_s1$r, cov, alp, pi_marg, brks.list[[1]])
      
      res <- -sum(obs.logL)
      
      print(res)
      
      return(res)
    }
    
    dim_phi1 <- R+1+len.eta1+len.beta
    res.s1 <- optim(par=rep(-1, dim_phi1), fn=obj.f,control = list(maxit= 1000),  method="L-BFGS-B", hessian=TRUE)
  }
  
  if(marg.type == "type-specific x nocovar"){
    
    J = max(indata_s1$type)
    indata_s1_j.list = lapply(1:J, function(j){
      subset(indata_s1, type==j)
    })
    
    ncol_Gmat <- dim(indata_s1[,-(1:7)])[2] 
    R = length(brks.list[[1]]) + 1
    
    obj.f <- function(phi1){
      
      #eta_1 = phi1[J*(R+1)+1:ncol_Gmat]
      
      res1 <- lapply(indata_s1_j.list, function(dt){
        
        XGmat <- matrix( as.numeric(as.matrix(dt[,-(1:7)])), ncol=ncol_Gmat, byrow=0)
        
        j = unique(dt$type)
        alp =  exp(phi1[(j-1)*(R+1)+(1:R)]); eta_0 = phi1[j*(R+1)]
        #pi_marg = expit.f(eta_0+XGmat %*% eta_1)
        pi_marg = expit.f(eta_0)
        
        beta = phi1[j*(R+1)+(1:length(adjust.cov.T))]
        cov = XGmat[,adjust.cov.T,drop=F] %*% beta
        
        obs.logL <- eval.obs.logL_adjusted_j.f(dt$l, dt$r, cov, alp, pi_marg, brks.list[[j]])
        res1 = ifelse(is.finite(obs.logL), obs.logL, 0)
        res0 = -sum(res1)
        #print(res0)
        res0 <- ifelse(is.finite(res0), res0, .Machine$double.xmax)
        return(res0)
      })
      
      res <- do.call("sum", res1)
      
      return(res)
    }
    
    dim_phi1 <- J*(R + 1)+length(adjust.cov.T)
    system.time(
      #res.s1 <- nlm(f=obj.f, p=rep(-1, dim_phi1),gradtol=1e-04, steptol=1e-05, hessian=TRUE)
      res.s1 <- optim(par=rep(-1, dim_phi1), fn=obj.f, control = list(maxit= 1000), method="L-BFGS-B", hessian = TRUE)
    )
  }
  
  #if(marg.type == "type-specific x covar-common"){
    
  #  J = max(indata_s1$type)
  #  indata_s1_j.list = lapply(1:J, function(j){
  #    subset(indata_s1, type==j)
  #  })
    
  #  ncol_Gmat <- dim(indata_s1[,-(1:7)])[2] 
  #  R = length(brks.list[[1]]) + 1
    
  #  obj.f <- function(phi1){
  #    
  #    eta_1 = phi1[J*(R+1)+1:length(adjust.cov.Z)]
      
  #    #print(eta_1)
      
  #    res1 <- lapply(indata_s1_j.list, function(dt){
        
  #      XGmat <- matrix( as.numeric(as.matrix(dt[,-(1:7)])), ncol=ncol_Gmat, byrow=0)
        
  #      j = unique(dt$type)
  #      alp =  exp(phi1[(j-1)*(R+1)+(1:R)]); eta_0 = phi1[j*(R+1)]
        
  #      if(is.null(adjust.cov.T)){
  #        cov=rep(0.0, dim(XGmat)[1])
  #      }else{
  #        beta = phi1[j*(R+1)+length(adjust.cov.Z)+(1:length(adjust.cov.T))]
  #        cov = XGmat[,adjust.cov.T,drop=F] %*% beta
  #      }
        
  #      pi_marg = expit.f(eta_0+XGmat[,adjust.cov.Z,drop=F] %*% eta_1)
        
  #      obs.logL <- eval.obs.logL_adjusted_j.f(dt$l, dt$r, cov, alp, pi_marg, brks.list[[j]])
  #      res2 <- ifelse(is.finite(obs.logL), obs.logL, 0)
  #      res0 = -sum(res2)
  #      res0
  #    })
      
  #    res <- do.call("sum", res1)
      
  #    print(phi1)
  #    print(res)
      
  #    return(res)
  #  }
    
  #  dim_phi1 <- J*(R + 1) + length(adjust.cov.Z)+length(adjust.cov.T)*(!is.null(adjust.cov.T))
  #  system.time(
  #    #res.s1 <- nlm(f=obj.f, p=rep(-1, dim_phi1),gradtol=1e-04, steptol=1e-05, hessian=TRUE)
  #    res.s1 <- optim(par=rep(-0.1, dim_phi1), fn=obj.f, 
  #                    control = list(maxit= 1000), method="L-BFGS-B", hessian = TRUE)
  #  )
    
  #}
  
  
  if(marg.type == "type-specific x covar-common"){
    
    obj.f <- function(phi1){
      
      #eta_1 = phi1[sum(R_vec+1)+1:len.eta1]
      eta_1 = phi1[J*(R+1)+1:len.eta1]
      
      if(is.null(adjust.cov.T)){
        beta=0
      }else{
        beta =  phi1[sum(R_vec+1)+len.eta1+(1:len.beta)]
      }
      
      
      res1 <- lapply(indata_s1.list, function(dt){
        
        j = unique(dt$type)
        alp =  exp(phi1[(j-1)*(R+1)+(1:R)]); eta_0 = phi1[j*(R+1)]
        #xgmat <- XGmat.list[[j]]
        
        XGmat <- matrix( as.numeric(as.matrix(dt[,-(1:7)])), ncol=ncol_Gmat, byrow=0)
        
        pi_marg = expit.f(eta_0+ XGmat[,adjust.cov.Z, drop=F] %*% eta_1)
        if(is.null(adjust.cov.T)){
          cov = XGmat[,1,drop=F] %*% beta
        }else{
          cov = XGmat[,adjust.cov.T,drop=F] %*% beta
        }
        
        obs.logL <- eval.obs.logL_adjusted_j.f(dt$l, dt$r,cov, alp, pi_marg, brks.list[[j]])
        #res2 <- ifelse(is.finite(obs.logL), obs.logL, 0)
        #res0 = -sum(res2)
        res0 = -sum(obs.logL)
        #print(res0)
        res0
      })
      
      res <- do.call("sum", res1)
      
      res <- ifelse(is.finite(res), res, .Machine$double.xmax)
      
      
      print(phi1); print(res)
      
      return(res)
    }
    
    dim_phi1 = J*(R + 1) +len.eta1+len.beta #sum(R_vec+1)+len.eta1+len.beta
    system.time(
      #res.s1 <- nlm(f=obj.f, p=rep(-1, dim_phi1),gradtol=1e-04, steptol=1e-05, hessian=TRUE)
      res.s1 <- optim(par= c(rep(log(0.01),2), -0.01, rep(log(0.01),2), -rep(0.01, 23)), #rep(-0.1, dim_phi1), 
                      fn=obj.f, 
                      control = list(maxit= 1000), method="L-BFGS-B", hessian = TRUE)
    )
    
  }
  
  score.s1 <- realdt_score_adjusted_phi1.f(gradstep=1e-06, phi1=res.s1$par, indata_s1, adjust.cov.Z, adjust.cov.T, brks.list, marg.type) 
  #hess.s1 <- res.s1$hessian
  #H.s1 <- (t(score.s1) %*% score.s1)
  #ase.s1 <- sqrt(diag(ginv( hess.s1 %*% ginv(H.s1) %*% hess.s1)))
  
  
  return(list(optim.res=res.s1,
              score=score.s1))
}


realdt_score_adjusted_phi1.f <- function(gradstep, phi1, indata_s1, adjust.cov.Z, adjust.cov.T,
                                brks.list, marg.type="common"){
  
  lik_phi1 <- realdt_lik_adjusted_phi1.f(phi1, indata_s1, adjust.cov.Z, adjust.cov.T, brks.list, marg.type)
  
  p1 = length(phi1)
  #score_mat <- matrix(NA, nrow=nsample, ncol=p1)
  
  score.list <- lapply(1:p1, function(x){
    
    phi1p = phi1
    phi1p[x] = phi1[x] + gradstep
    
    profile_x <- realdt_lik_adjusted_phi1.f(phi1p, indata_s1, adjust.cov.Z, adjust.cov.T, brks.list, marg.type)
    tmp <- (profile_x - lik_phi1)/gradstep
    
    tmp
  })
  
  score_mat <- do.call("cbind", score.list)
  
  return(score_mat)
}


realdt_hess_adjusted_phi1.f <- function(gradstep, phi1, indata_s1, adjust.cov.Z, adjust.cov.T,
                                        brks.list, marg.type="common"){
  
  p = length(phi1)
  
  score_vec_0 <- colSums(realdt_score_adjusted_phi1.f(gradstep, phi1, indata_s1, adjust.cov.Z, adjust.cov.T, brks.list, marg.type))
  
  sv_p.list <- lapply(1:p, function(x){
    
    vartheta_p = phi1 
    vartheta_p[x] = phi1[x] + gradstep
    
    score_vec_p <- colSums(realdt_score_adjusted_phi1.f(gradstep, vartheta_p, indata_s1, adjust.cov.Z, adjust.cov.T, brks.list, marg.type))
    score_vec_p
  } )
  
  Hess <- matrix(0.0, p, p)
  for(j in 1:p){
    
    Hess[j,] <- (sv_p.list[[j]] - score_vec_0)/gradstep
    
  }
  
  
  return(Hess)
  
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
      eta_11 = phi1[c(R1+R2+1+(1:len.eta1))]
      eta_21 = phi1[c(R1+R2+2+len.eta1+(1:len.eta1))]
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
  }else if(marg.type=="type-specific x adj.T.cov"){
    
    alp_1 = exp(phi1[1:R1]);  
    eta_10 = phi1[R1+1];   
    #eta_11 = phi1[c(R1+1+(1:ncol_Gmat))]; 
    alp_2 = exp(phi1[R1+1+(1:R2)]); 
    eta_20 = phi1[R1+R2+2];
    
    if(len.eta1>0){
      eta_21 = eta_11 = phi1[c(R1+R2+2+(1:len.eta1))]
      pi_marg_1 = expit.f(eta_10+XGmat_j1[,adjust.cov.Z, drop=F] %*% eta_11)
      pi_marg_2 = expit.f(eta_20+XGmat_j2[,adjust.cov.Z, drop=F] %*% eta_21)
    }else{
      pi_marg_1 = expit.f(rep(eta_10, dim(XGmat_j1)[1]))
      pi_marg_2 = expit.f(rep(eta_20, dim(XGmat_j2)[1]))
    }
    
    if(len.beta>0){
      beta_11 = phi1[c(R1+R2+2+(len.eta1)+(1:len.beta))]
      beta_21 = phi1[c(R1+R2+2+(len.eta1)+len.beta+(1:len.beta))]
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

pseudoScoreTest_adj.f <- function(gradstep, chi.df, phi1.expand, indata_s1, XGmat, adjust.cov.Z, adjust.cov.T,
                              brks.list, marg.type="type-specific x covar-common"){
  
  len.beta <- ifelse(is.null(adjust.cov.T), 0, length(adjust.cov.T))
  
  s2 <- score_expand_phi1.f(gradstep=1e-06, phi1=phi1.expand, indata_s1, XGmat, 
                            adjust.cov.Z, adjust.cov.T, brks.list, marg.type) 
  nsample = length(unique(indata_s1$id))
  J2 <- (t(s2) %*% s2)/nsample
  
  est.pst <- (colSums(s2) %*% ginv(J2) %*% colSums(s2))/nsample
  
  pval=1-pchisq(est.pst, df=chi.df)
  
  res.list <- list(pst.stat=est.pst,
                   pval=pval)
  return(res.list)
}


#XGmat <- dt_s1[,c(paste0("x", adjust.cov.Z))]
#pseudoScoreTest_adj.f(1e-06, 11, c(res.s1$par, res.s1$par[(length(res.s1$par)-len.beta+1):length(res.s1$par)]), 
#                      indata_s1, XGmat, adjust.cov.Z, adjust.cov.T, 
#                      brks.list, marg.type="type-specific x adj.T.cov")





