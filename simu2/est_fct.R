
# 1. working-indep composite lik.
# 2. pairwise composite lik.

# bivariate normal: CDF
my_pbivnorm <- function(x, y, rho){
  
  x.Inf <- is.infinite(x)
  y.Inf <- is.infinite(y)
  
  res0 <- pbivnorm(x, y, rho)
  pos.Inf.ind <- (x.Inf + y.Inf)>=1
  pos.Inf <- which((pos.Inf.ind))
  
  if(any(pos.Inf.ind)){
    
    x1 <- x[pos.Inf]
    y1 <- y[pos.Inf]
    
    res1 <- pnorm(x1)*pnorm(y1)
    
    res0[pos.Inf] <- res1
    
  }
  
  return( res0)
  
}

# working-indep loglik
Lcomp_1j <- function(l, r, alp, brks){
  #res <- ppch(r, cuts=brks, levels=alp) - ppch(l, cuts=brks, levels=alp)
  res <- vppwc(r, cuts=brks, levels=alp, 1, 0) - vppwc(l, cuts=brks, levels=alp, 1, 0)
  return(as.double(res))
}

eval.obs.logL_j.f <- function(l, r, alp, pi_marg, brks){
  
  comp.L_1j <- Lcomp_1j(l, r, alp, brks)
  use.finite <- is.finite(r)
  pi_marg_0 <- 1- pi_marg
  if(any(use.finite)){ pi_marg_0[use.finite] <- 0 }
  
  res <- log(comp.L_1j*pi_marg + pi_marg_0)
  
  res <- ifelse(is.nan(res), .Machine$double.xmin, res)
  return(res)
}

# --------- 1. estimate stage I parameters phi1 = c(log(alp_vec), eta_vec)
# --- 1.1. working-indep lik evaluation using Rcpp functions
# ------- evaluation of loglik is using R instead of Rcpp as the former is much more faster
estParam_s1.f <- function(indata_s1.list, brks.list, marg.type=NULL){
  
  #if(!(marg.type %in% c("common", "type-specific", "gene-common"))){
    
  #  stop("marg.type has to be common, type-specific, or gene-common !")
    
  #}
  
  indata_s1 = do.call("rbind", indata_s1.list)
  
  
  if(is.null(marg.type)){
    marg.type = "type-specific"
  }
  
  if(marg.type == "common"){
    
    R <- length(brks.list[[1]])+1
    
    
    obj.f <- function(phi1){
      
      alp = exp(phi1[1:R]);
      eta = phi1[-(1:R)]; 
      pi_marg = expit.f(eta[1]+eta[2]*indata_s1$x+eta[3]*indata_s1$g)
      
      obs.logL <- eval.obs.logL_j.f(indata_s1$l, indata_s1$r, alp, pi_marg, brks.list[[1]])
      
      res <- -sum(obs.logL)
      
      return(res)
    }
    
    dim_phi1 <- length(brks_vec)+4
    res.s1 <- nlm(f=obj.f, p=rep(0.05,dim_phi1),gradtol=1e-06, steptol=1e-06)
    
  }
  
  if(marg.type == "type-specific"){
    
    indata_s1_j1 = subset(indata_s1, type==1)
    indata_s1_j2 = subset(indata_s1, type==2)
    indata_s1_j3 = subset(indata_s1, type==3)
    
    R1 = length(brks.list[[1]]) + 1; R2 = length(brks.list[[2]]) + 1; R3 = length(brks.list[[3]]) + 1
    R = R1 + R2 + R3
    
    obj.f <- function(phi1){
      
      alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[R1+1+(1:R2)]); alp_3 = exp(phi1[R1+R2+2+(1:R3)])
      eta_10 = phi1[R1+1]; eta_20 = phi1[R1+R2+2];  eta_30 = phi1[R+3]
      eta_x = phi1[R+4];
      eta_g_1 = phi1[R+5]; eta_g_2 = phi1[R+6]; eta_g_1 = phi1[R+7]
      pi_marg_1 = expit.f(eta_10+eta_x*indata_s1_j1$x+eta_g_1*indata_s1_j1$g)
      pi_marg_2 = expit.f(eta_20+eta_x*indata_s1_j2$x+eta_g_2*indata_s1_j2$g)
      pi_marg_3 = expit.f(eta_30+eta_x*indata_s1_j3$x+eta_g_3*indata_s1_j3$g)
      
      obs.logL_1 <- eval.obs.logL_j.f(indata_s1_j1$l, indata_s1_j1$r, alp_1, pi_marg_1, brks.list[[1]])
      obs.logL_2 <- eval.obs.logL_j.f(indata_s1_j2$l, indata_s1_j2$r, alp_2, pi_marg_2, brks.list[[2]])
      obs.logL_3 <- eval.obs.logL_j.f(indata_s1_j3$l, indata_s1_j3$r, alp_3, pi_marg_3, brks.list[[2]])
      
      res <- -(sum(obs.logL_1)+sum(obs.logL_2)+sum(obs.logL_3))
      
      return(res)
    }
    
    dim_phi1 <- R + 3 + 1 + 3
    res.s1 <- nlm(f=obj.f, p=rep(0.05,dim_phi1),gradtol=1e-06, steptol=1e-06, hessian=TRUE)
    
  }
  
  if(marg.type == "gene-common"){
    
    indata_s1_j1 = subset(indata_s1, type==1)
    indata_s1_j2 = subset(indata_s1, type==2)
    indata_s1_j3 = subset(indata_s1, type==3)
     
    R1 = length(brks.list[[1]]) + 1; R2 = length(brks.list[[2]]) + 1; R3 = length(brks.list[[3]]) + 1
    R = R1 + R2 + R3
    
    obj.f <- function(phi1){
      
      alp_1 = exp(phi1[1:R1]); alp_2 = exp(phi1[R1+1+(1:R2)]); alp_3 = exp(phi1[R1+R2+2+(1:R3)])
      eta_10 = phi1[R1+1]; eta_20 = phi1[R1+R2+2];  eta_30 = phi1[R+3]
      eta_xg = phi1[R+3+c(1,2)]
      pi_marg_1 = expit.f(eta_10+eta_xg[1]*indata_s1_j1$x+eta_xg[2]*indata_s1_j1$g)
      pi_marg_2 = expit.f(eta_20+eta_xg[1]*indata_s1_j2$x+eta_xg[2]*indata_s1_j2$g)
      pi_marg_3 = expit.f(eta_30+eta_xg[1]*indata_s1_j3$x+eta_xg[2]*indata_s1_j3$g)
      
      obs.logL_1 <- eval.obs.logL_j.f(indata_s1_j1$l, indata_s1_j1$r, alp_1, pi_marg_1, brks.list[[1]])
      obs.logL_2 <- eval.obs.logL_j.f(indata_s1_j2$l, indata_s1_j2$r, alp_2, pi_marg_2, brks.list[[2]])
      obs.logL_3 <- eval.obs.logL_j.f(indata_s1_j3$l, indata_s1_j3$r, alp_3, pi_marg_3, brks.list[[2]])
      
      res <- -(sum(obs.logL_1)+sum(obs.logL_2)+sum(obs.logL_3))
      
      return(res)
    }
    
    dim_phi1 <- R + 3 + 2 
    system.time(
    res.s1 <- nlm(f=obj.f, p=rep(0.01, dim_phi1),gradtol=1e-06, steptol=1e-06, hessian=TRUE)
    )
    
  }
  
  return(res.s1)
}


# --------- 2. estimate stage II parameters phi2 = c(log(or_pair), atan(cor_pair)*2/pi)
# --- 2.1 
# ------- (1) need to prepare phi1_est.list (length of J) first
estParam_s2.f <- function(phi1_est.list, brks.list, 
                          indata_s2.list, pair.type=NULL){
  
  if(is.null(pair.type)){
    pair.type = "jjp"
  }
  
  if(pair.type=="jjp"){
    
    res_s2.list <- mclapply(indata_s2.list, function(dt){
      
      j = unique(dt$type_j)
      jp = unique(dt$type_jp)
      
      brks_j = brks.list[[j]]
      brks_jp = brks.list[[jp]]
      
      phi1_j = phi1_est.list[[j]]
      phi1_jp = phi1_est.list[[jp]]
      
      
      obj.f <- function(phi2){
        
        #newlogCLPair_jkjpkp()
        
        #res <- newlogCLPair_jkjpkp(phi1_j, phi1_jp, phi2, dt$l_j, dt$r_j,
        #                        dt$l_jp, dt$r_jp, as.matrix(dt[, c("x","g")]), 
        #                        brks_j, brks_jp)
        
        indata_dt <- get_indata_s2_jjp(phi1_j, phi1_jp, dt$l_j, dt$r_j,
                                       dt$l_jp, dt$r_jp, 
                                       as.matrix(dt[, c("x","g")]), as.matrix(dt[, c("x","g")]),
                                       brks_j, brks_jp)
        
        sigma_jjp = atan(phi2[2])*2/pi
        gamma = phi2[1]
        
        #t1 = Sys.time()
        #d2_jjp <- bi.mvn.me(cbind(indata_dt$invPhisurvlt_j, indata_dt$invPhisurvlt_jp), sigma_jjp) - 
        #  bi.mvn.me(cbind(indata_dt$invPhisurvlt_j, indata_dt$invPhisurvrt_jp), sigma_jjp) - 
        # bi.mvn.me(cbind(indata_dt$invPhisurvrt_j, indata_dt$invPhisurvlt_jp), sigma_jjp) +
        #  bi.mvn.me(cbind(indata_dt$invPhisurvrt_j, indata_dt$invPhisurvrt_jp), sigma_jjp)
        d2_jjp <- (my_pbivnorm(indata_dt$invPhisurvlt_j, indata_dt$invPhisurvlt_jp, sigma_jjp) - 
                     my_pbivnorm(indata_dt$invPhisurvlt_j, indata_dt$invPhisurvrt_jp, sigma_jjp) - 
                     my_pbivnorm(indata_dt$invPhisurvrt_j, indata_dt$invPhisurvlt_jp, sigma_jjp) +
                     my_pbivnorm(indata_dt$invPhisurvrt_j, indata_dt$invPhisurvrt_jp, sigma_jjp) )
        #t2 = Sys.time()
        #print(t2-t1)
        
        #t1 = Sys.time()
        res <- newevalLogCLPair_jjp(gamma, d2_jjp, is.finite(dt$r_j), is.finite(dt$r_jp), 
                                    indata_dt$pi_marg_j, indata_dt$pi_marg_jp, 
                                    indata_dt$survlt_j, indata_dt$survrt_j,
                                    indata_dt$survlt_jp, indata_dt$survrt_jp)
        #t2 = Sys.time()
        #print(t2-t1)
        
        #print(sum(res))
        
        return(-sum(res))
        
      }
      
      #t1 = Sys.time()
      system.time(
      res.s2 <- nlm(obj.f, p=rep(0.01,2), gradtol=1e-06, steptol=1e-06, hessian=TRUE)
      )
      #t2 = Sys.time()
      #print(t2-t1)
      #res.s2 <- optim(par=rep(0.05, 2), fn=obj.f, method="BFGS")
      res.s2
    }, mc.cores=6)
    
  }
  
  return(res_s2.list)
  
}


