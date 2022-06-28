
# Numerical calculations for inference components
# score functions: working-indep. score, pairwise score
# sub Hessian matrix for the partial(phi2)partial(phi1)

lik_phi1.f <- function(phi1, indata_s1, brks.list, marg.type="gene-common"){
  
  indata_s1_j1 = subset(indata_s1, type==1)
  indata_s1_j2 = subset(indata_s1, type==2)
  indata_s1_j3 = subset(indata_s1, type==3)
  
  R1 = length(brks.list[[1]]) +1; R2 = length(brks.list[[2]]) +1; R3 = length(brks.list[[3]]) +1
  R = R1 + R2 + R3
  
  if(marg.type=="gene-common"){
    
    alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[R1+1+(1:R2)]); alp_3 = exp(phi1[R1+R2+2+(1:R3)])
    eta_10 = phi1[R1+1]; eta_20 = phi1[R1+R2+2];  eta_30 = phi1[R+3]
    eta_xg = phi1[R+3+c(1,2)]
    pi_marg_1 = expit.f(eta_10+eta_xg[1]*indata_s1_j1$x+eta_xg[2]*indata_s1_j1$g)
    pi_marg_2 = expit.f(eta_20+eta_xg[1]*indata_s1_j2$x+eta_xg[2]*indata_s1_j2$g)
    pi_marg_3 = expit.f(eta_30+eta_xg[1]*indata_s1_j3$x+eta_xg[2]*indata_s1_j3$g)
    
    obs.logL_1 <- eval.obs.logL_j.f(indata_s1_j1$l, indata_s1_j1$r, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_j.f(indata_s1_j2$l, indata_s1_j2$r, alp_2, pi_marg_2, brks.list[[2]])
    obs.logL_3 <- eval.obs.logL_j.f(indata_s1_j3$l, indata_s1_j3$r, alp_3, pi_marg_3, brks.list[[3]])
    
    res <- c(obs.logL_1, obs.logL_2, obs.logL_3)
    
  }else if(marg.type=="type-specific"){
    
    alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[R1+1+(1:R2)]); alp_3 = exp(phi1[R1+R2+2+(1:R3)])
    eta_10 = phi1[R1+1]; eta_20 = phi1[R1+R2+2];  eta_30 = phi1[R+3]
    eta_x = phi1[R+4]; eta_g_1 = phi1[R+5]; eta_g_2 = phi1[R+6]; eta_g_3 = phi1[R+7]
    
    pi_marg_1 = expit.f(eta_10+eta_x*indata_s1_j1$x+eta_g_1*indata_s1_j1$g)
    pi_marg_2 = expit.f(eta_20+eta_x*indata_s1_j2$x+eta_g_2*indata_s1_j2$g)
    pi_marg_3 = expit.f(eta_30+eta_x*indata_s1_j3$x+eta_g_3*indata_s1_j3$g)
    
    obs.logL_1 <- eval.obs.logL_j.f(indata_s1_j1$l, indata_s1_j1$r, alp_1, pi_marg_1, brks.list[[1]])
    obs.logL_2 <- eval.obs.logL_j.f(indata_s1_j2$l, indata_s1_j2$r, alp_2, pi_marg_2, brks.list[[2]])
    obs.logL_3 <- eval.obs.logL_j.f(indata_s1_j3$l, indata_s1_j3$r, alp_3, pi_marg_3, brks.list[[3]])
    
    res <- c(obs.logL_1, obs.logL_2, obs.logL_3)
    
  }

  
  re_indata_s1 <- as.data.frame(rbind(indata_s1_j1, indata_s1_j2, indata_s1_j3))
  
  res <- aggregate(res, list(re_indata_s1$id), sum, drop = T, simplify = T)[, -1]
  
  return(res)
  
}


score_phi1.f <- function(gradstep, phi1, nsample, indata_s1,
                         brks.list, marg.type="gene-common"){
  
  lik_phi1 <- lik_phi1.f(phi1, indata_s1, brks.list, marg.type)
  
  p1 = length(phi1)
  #score_mat <- matrix(NA, nrow=nsample, ncol=p1)
  
  score.list <- lapply(1:p1, function(x){
    
    phi1p = phi1
    phi1p[x] = phi1[x] + gradstep
    
    profile_x <- lik_phi1.f(phi1p, indata_s1, brks.list, marg.type=marg.type)
    tmp <- (profile_x - lik_phi1)/gradstep
    
    tmp
  })
  
  score_mat <- do.call("cbind", score.list)
  
  return(score_mat)
}

lik_phi2_jjp.f <- function(phi2, df_s2_jjp, phi1_j, phi1_jp, brks_j, brks_jp){
  
  indata_dt <- get_indata_s2_jjp(phi1_j, phi1_jp, df_s2_jjp$l_j, df_s2_jjp$r_j,
                                 df_s2_jjp$l_jp, df_s2_jjp$r_jp, 
                                 as.matrix(df_s2_jjp[, c("x","g")]), 
                                 as.matrix(df_s2_jjp[, c("x","g")]),
                                 brks_j, brks_jp)
  
  sigma_jjp = atan(phi2[2])*2/pi
  gamma = phi2[1]
  
  d2_jjp <- (my_pbivnorm(indata_dt$invPhisurvlt_j, indata_dt$invPhisurvlt_jp, sigma_jjp) - 
               my_pbivnorm(indata_dt$invPhisurvlt_j, indata_dt$invPhisurvrt_jp, sigma_jjp) - 
               my_pbivnorm(indata_dt$invPhisurvrt_j, indata_dt$invPhisurvlt_jp, sigma_jjp) +
               my_pbivnorm(indata_dt$invPhisurvrt_j, indata_dt$invPhisurvrt_jp, sigma_jjp) )

  res <- newevalLogCLPair_jjp(gamma, d2_jjp, is.finite(df_s2_jjp$r_j), is.finite(df_s2_jjp$r_jp), 
                              indata_dt$pi_marg_j, indata_dt$pi_marg_jp, 
                              indata_dt$survlt_j, indata_dt$survrt_j,
                              indata_dt$survlt_jp, indata_dt$survrt_jp)
  
  res <- aggregate(res, list(df_s2_jjp$id), sum, drop = T, simplify = T)[, -1]
  
  return(res)
  
}

score_phi2_jjp.f <- function(gradstep, phi2_jjp, nsample, df_s2_jjp,
                             phi1_j, phi1_jp, brks_j, brks_jp){
  
  lik_phi2_jjp <- lik_phi2_jjp.f(phi2_jjp, df_s2_jjp, phi1_j, phi1_jp, brks_j, brks_jp)
  
  score.list <- lapply(1:2, function(x){
    
    phi2p = phi2_jjp
    phi2p[x] = phi2_jjp[x] + gradstep
    
    profile_x <- lik_phi2_jjp.f(phi2p, df_s2_jjp, phi1_j, phi1_jp, brks_j, brks_jp)
    
    tmp <- (profile_x - lik_phi2_jjp)/gradstep
    
    tmp
  })
  
  score_mat <- do.call("cbind", score.list)
  
  return(score_mat)
  
}


Hess_IP_jjp.f <- function(gradstep, phi2_jjp, nsample, df_s2_jjp,
                         j, jp, phi1, indata_s1, brks.list, marg.type){
  
  R1 = length(brks.list[[1]]) +1; R2 = length(brks.list[[2]]) +1; R3 = length(brks.list[[3]]) +1
  R = R1 + R2 + R3
  
  if(marg.type == "gene-common"){
    phi1_1 = phi1[c(1:(R1+1), R+3+(1:2))]
    phi1_2 = phi1[c(R1+1+(1:(R2+1)), R+3+c(1,2))]
    phi1_3 = phi1[c(R1+R2+2+(1:(R3+1)), R+3+c(1,2))]
  }else if(marg.type == "type-specific"){
    
    phi1_1 = phi1[c(1:(R1+1), R+3+c(1,2))]
    phi1_2 = phi1[c(R1+1+(1:(R2+1)), R+3+c(1,3))]
    phi1_3 = phi1[c(R1+R2+2+(1:(R3+1)), R+3+c(1,4))]
    
  }
  
  
  phi1.list <- list(phi1_1, phi1_2, phi1_3)
  
  phi1_j = phi1.list[[j]]; phi1_jp = phi1.list[[jp]]
  brks_j = brks.list[[j]]; brks_jp = brks.list[[jp]]
  
  #lik_phi1 <- sum(lik_phi1.f(phi1, indata_s1, R1, R2, R3, marg.type))
  
  lik_phi2_jjp <- sum(lik_phi2_jjp.f(phi2_jjp, df_s2_jjp, phi1_j, phi1_jp, brks_j, brks_jp))
  
  p1 = length(phi1)
  #score_mat <- matrix(NA, nrow=nsample, ncol=p1)
  
  pf_phi1.list <- lapply(1:p1, function(x){
    
    phi1p = phi1
    phi1p[x] = phi1[x] + gradstep
    
    if(marg.type=="gene-common"){
      phi1p_1 = phi1p[c(1:(R1+1), R+3+(1:2))]
      phi1p_2 = phi1p[c(R1+1+(1:(R2+1)), R+3+c(1,2))]
      phi1p_3 = phi1p[c(R1+R2+2+(1:(R3+1)), R+3+c(1,2))]
    }else if(marg.type == "type-specific"){
      
      phi1p_1 = phi1p[c(1:(R1+1), R+3+(1:2))]
      phi1p_2 = phi1p[c(R1+1+(1:(R2+1)), R+3+c(1,3))]
      phi1p_3 = phi1p[c(R1+R2+2+(1:(R3+1)), R+3+c(1,4))]
    }
    
    
    
    phi1p_est.list <- list(phi1p_1, phi1p_2, phi1p_3)
    
    phi1p_j = phi1p_est.list[[j]]; phi1p_jp = phi1p_est.list[[jp]]
    
    tmp <- sum(lik_phi2_jjp.f(phi2_jjp, df_s2_jjp, phi1p_j, phi1p_jp, brks_j, brks_jp))
    #tmp <- sum(lik_phi1.f(phi1p, indata_s1, R1, R2, R3, marg.type="gene-common"))
    
    tmp
  })
  
  pf_phi2_jjp.list <- lapply(1:2, function(x){
    
    phi2p = phi2_jjp
    phi2p[x] = phi2_jjp[x] + gradstep
    
    tmp <- sum(lik_phi2_jjp.f(phi2p, df_s2_jjp, phi1_j, phi1_jp, brks_j, brks_jp))
    tmp
  })
  
  pf_phi1_vec <- as.vector(unlist(pf_phi1.list))
  pf_phi2_vec <- as.vector(unlist(pf_phi2_jjp.list))
  
  Hess <- matrix(0.0, 2, p1)
  for(x in 1:2){
    for(y in 1:p1){
      
      phi2p = phi2_jjp ; phi1p = phi1
      phi2p[x] = phi2_jjp[x] + gradstep; phi1p[y] = phi1[y] + gradstep
      
      if(marg.type=="gene-common"){
        phi1p_1 = phi1p[c(1:(R1+1), R+3+(1:2))]
        phi1p_2 = phi1p[c(R1+1+(1:(R2+1)), R+3+c(1,2))]
        phi1p_3 = phi1p[c(R1+R2+2+(1:(R3+1)), R+3+c(1,2))]
      }else if(marg.type == "type-specific"){
        
        phi1p_1 = phi1p[c(1:(R1+1), R+3+(1:2))]
        phi1p_2 = phi1p[c(R1+1+(1:(R2+1)), R+3+c(1,3))]
        phi1p_3 = phi1p[c(R1+R2+2+(1:(R3+1)), R+3+c(1,4))]
      }
      
      phi1p_est.list <- list(phi1p_1, phi1p_2, phi1p_3)
      
      phi1p_j = phi1p_est.list[[j]]; phi1p_jp = phi1p_est.list[[jp]]
      
      tmp <- sum(lik_phi2_jjp.f(phi2p, df_s2_jjp, phi1p_j, phi1p_jp, brks_j, brks_jp))
      
      Hess[x,y] <- (tmp- pf_phi1_vec[y] + lik_phi2_jjp - pf_phi2_vec[x])/(gradstep*gradstep) 
      
    }
  }
  
  return(Hess)
}


ts_asvar.f <- function(Hess_s1, Hess_phi2.list, 
                       gradstep, nsample, 
                       phi1, brks.list, indata_s1, 
                       phi2.list, df_s2.list,
                       marg.type, num.cores=2){
  
  p1 = length(phi1)
  
  jjp.mat <- rbind(c(1,1), c(2,2), c(3,3),
                   c(1,2), c(1,3), c(2,3))
  
  R1 = length(brks.list[[1]]) +1; R2 = length(brks.list[[2]]) +1; R3 = length(brks.list[[3]]) +1
  R = R1 + R2 + R3
  
  # score_phi1.mat
  score_phi1 <- score_phi1.f(gradstep, phi1, nsample, indata_s1, brks.list, marg.type)
  
  # score_phi2.mat
  
  if(marg.type == "gene-common"){
    phi1_1 = phi1[c(1:(R1+1), R+3+(1:2))]
    phi1_2 = phi1[c(R1+1+(1:(R2+1)), R+3+c(1,2))]
    phi1_3 = phi1[c(R1+R2+2+(1:(R3+1)), R+3+c(1,2))]
  }else if(marg.type == "type-specific"){
    
    phi1_1 = phi1[c(1:(R1+1), R+3+c(1,2))]
    phi1_2 = phi1[c(R1+1+(1:(R2+1)), R+3+c(1,3))]
    phi1_3 = phi1[c(R1+R2+2+(1:(R3+1)), R+3+c(1,4))]
    
  }
  phi1.list <- list(phi1_1, phi1_2, phi1_3)
  
  score_phi2.list <- mclapply(1:6, function(x){
    
    j = jjp.mat[x,1]; jp = jjp.mat[x,2]
    score_phi2_jjp.f(gradstep, phi2.list[[x]],
                     nsample, df_s2.list[[x]], phi1.list[[j]],
                     phi1.list[[jp]], brks.list[[j]], brks.list[[jp]])
  }, mc.cores=num.cores)
  
  score_phi2 <- do.call("cbind", score_phi2.list)
  
  # Hess_IP: 2 x p1
  Hess_IP.list <- mclapply(1:6, function(x){
    
    j = jjp.mat[x,1]; jp = jjp.mat[x,2]
    Hess_IP_jjp.f(gradstep, phi2.list[[x]],
                  nsample, df_s2.list[[x]], j, jp, 
                  phi1, indata_s1, brks.list, marg.type)
  }, mc.cores=num.cores)
  
  Hess_IP <- do.call("rbind", Hess_IP.list)
  
  # Hess_jj
  #for(x in 1:6){
  #  Hess_s2[x*2+(c(-1,0)), x*2+(c(-1,0))] <- Hess_phi2.list[[x]]
  #}
  
  Hess_s2.list <- lapply(1:6, function(x){
    Hess_s2 <- matrix(0, nrow=2, ncol=12)
    Hess_s2[, x*2+c(-1,0)] <- Hess_phi2.list[[x]]
    Hess_s2
  })
  Hess_s2 <- do.call("rbind", Hess_s2.list)
  
  # A
  A11 = -Hess_s1 
  A21 = -Hess_IP
  A22 = -Hess_s2
  #A = rbind(A1, A2)
  score_phi1 <- as.matrix(score_phi1); score_phi2 = as.matrix(score_phi2)
  #score_comb <- cbind(score_phi1, score_phi2)
  #B = t(score_comb) %*% score_comb
  B11 = t(score_phi1) %*% score_phi1
  B22 = t(score_phi2) %*% score_phi2
  
  #ts_asvar = ginv(A) %*% B %*% t(ginv(A))
  
  ts_asvar_phi1 = ginv(A11) %*% B11 %*% t(ginv(A11))
  
  ts_asvar_phi2 = ginv(A22) %*% (B22 - A21%*% ginv(A11) %*% t(A21)) %*% t(ginv(A22))
  
  #sqrt(diag(ginv(B22)))
  #sqrt(diag(ginv(A22)))
  
  ts_ase_phi1 = sqrt(diag(ts_asvar_phi1))
  ts_ase_phi2 = sqrt(diag(ts_asvar_phi2))
  
  if(any(is.na(ts_ase_phi1))){
    print(phi1)
    stop("phi1.ase is NaN!")
  }
  
  if(any(is.na(ts_ase_phi2))){
    print(phi1)
    print(phi2.list)
    stop("phi2.ase is NaN!")
  }
  
  res <- list(phi1=ts_ase_phi1, phi2=ts_ase_phi2)
  
  return(res)
}




