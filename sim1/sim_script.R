
# ------------------------------------------------------------------------------------------
# script of the 1st set of simulation studies
# used to generate results presented in Table 2 and Web Table 1 
# ------------------------------------------------------------------------------------------

library(Rcpp)
library(Rsolnp)
library(MASS)
library(NlcOptim)
library(eha)
library(parallel)
library(mnormt)
library(poisson)
library(pbivnorm)
sourceCpp("commonf.cpp")
sourceCpp("loglik.cpp")
source("dataGen_fct.R")
source("est_fct.R")
source("ts_inference_fct.R")


# script function for simulation studies, include:
# 1. data generation
# 2. estimation and inference
# 3. record results in .dat files
J3_twoStage_script.f <- function(nsim, nsample, 
                                 prob.g, prob.x,
                                 K1, K2, K3, 
                                 pi1, pi2, pi3,
                                 eta_11_vec, eta_21_vec, eta_31_vec,
                                 SurvP.ctr1, SurvP.ctr2, SurvP.ctr3,
                                 #alp_11_vec, alp_21_vec, alp_31_vec,
                                 kappa_1, kappa_2, kappa_3,
                                 brks_1_vec, brks_2_vec, brks_3_vec,
                                 or_11, or_22, or_33, or_12, or_13, or_23, 
                                 tau_11, tau_22, tau_33, tau_12, tau_13, tau_23,
                                 visit.rate,
                                 marg.type = NULL, pair.type = NULL
                                 ){
  
  # -------- Parameter preparation --------- #
  # prepare stage I params: alp_j_vec, eta_j_vec
  eta_10 = sol.eta_0.f(pi1, eta_11_vec, c(prob.x, prob.g))
  eta_20 = sol.eta_0.f(pi2, eta_21_vec, c(prob.x, prob.g))
  eta_30 = sol.eta_0.f(pi3, eta_31_vec, c(prob.x, prob.g))
  
  #alp_10 = sol.alp_0.f(SurvP.ctr1, 1, brks_1_vec, kappa_1)
  #alp_20 = sol.alp_0.f(SurvP.ctr2, 1, brks_2_vec, kappa_2)
  #alp_30 = sol.alp_0.f(SurvP.ctr3, 1, brks_3_vec, kappa_3)
  
  eta_1_vec = c(eta_10, eta_11_vec); eta_2_vec = c(eta_20, eta_21_vec); eta_3_vec = c(eta_30, eta_31_vec)
  #alp_1_vec = alp_10*exp(kappa_1*(0:length(brks_1_vec)))
  #alp_2_vec = alp_20*exp(kappa_2*(0:length(brks_2_vec)))
  #alp_3_vec = alp_10*exp(kappa_3*(0:length(brks_3_vec)))
  
  scale_1 = sol.scale.f(SurvP.ctr1, 1, kappa_1)
  scale_2 = sol.scale.f(SurvP.ctr2, 1, kappa_2)
  scale_3 = sol.scale.f(SurvP.ctr3, 1, kappa_3)
  
  # prepare stage II params: pi_pair_mat, cor_mat, Pjoint_Z_g_XG
  pi_pair_mat.res = create_pi_pair_mat.f(K1, K2, K3, pi1, pi2, pi3, or_11, or_22, or_33, or_12, or_13, or_23)
  pi_pair_mat = pi_pair_mat.res$pi_pair_mat
  cor_mat = create_cor_mat.f(K1, K2, K3, tau_11, tau_22, tau_33, tau_12, tau_13, tau_23)
  
  pjoint_Z_g_11 = J3_Pjoint_Z_g_XG.f(K1, K2, K3, eta_1_vec, eta_2_vec, eta_3_vec, 1, 1, 
                                     pi_pair_mat,  
                                     pi_pair_mat.res$pi11, pi_pair_mat.res$pi22, pi_pair_mat.res$pi33, 
                                     pi_pair_mat.res$pi12, pi_pair_mat.res$pi13, pi_pair_mat.res$pi23)
  if(any(pjoint_Z_g_11<0)){stop("pjoint_Z_g_11 is not legistimate !")}
  
  pjoint_Z_g_10 = J3_Pjoint_Z_g_XG.f(K1, K2, K3, eta_1_vec, eta_2_vec, eta_3_vec, 1, 0, pi_pair_mat, 
                                     pi_pair_mat.res$pi11, pi_pair_mat.res$pi22, 
                                     pi_pair_mat.res$pi33, pi_pair_mat.res$pi12, 
                                     pi_pair_mat.res$pi13, pi_pair_mat.res$pi23)
  if(any(pjoint_Z_g_10<0)){stop("pjoint_Z_g_10 is not legistimate !")}
  
  pjoint_Z_g_01 = J3_Pjoint_Z_g_XG.f(K1, K2, K3, eta_1_vec, eta_2_vec, eta_3_vec, 0, 1, pi_pair_mat,
                                     pi_pair_mat.res$pi11, pi_pair_mat.res$pi22, 
                                     pi_pair_mat.res$pi33, pi_pair_mat.res$pi12, 
                                     pi_pair_mat.res$pi13, pi_pair_mat.res$pi23)
  if(any(pjoint_Z_g_01<0)){stop("pjoint_Z_g_01 is not legistimate !")}
  
  pjoint_Z_g_00 = J3_Pjoint_Z_g_XG.f(K1, K2, K3, eta_1_vec, eta_2_vec, eta_3_vec, 0, 0, pi_pair_mat,
                                     pi_pair_mat.res$pi11, pi_pair_mat.res$pi22, 
                                     pi_pair_mat.res$pi33, pi_pair_mat.res$pi12, 
                                     pi_pair_mat.res$pi13, pi_pair_mat.res$pi23)
  if(any(pjoint_Z_g_00<0)){stop("pjoint_Z_g_00 is not legistimate !")}
  
  
  K= K1+K2+K3
  l <- rep(list(0:1), K)
  Zmat <- expand.grid(l)
  reZmat = rbind(Zmat[rowSums(Zmat)<=1,],
                 Zmat[rowSums(Zmat)==2,],
                 Zmat[rowSums(Zmat)>=3,])
  
  R1 = length(brks_1_vec)+1; R2 = length(brks_2_vec)+1; R3 = length(brks_3_vec)+1
  
  # -------- Record simulation results --------- #
  
  file.tx <- paste("2stg.n", nsample,"J3", "K1", K1, "K2", K2, "K3", K3,  
                   "prob.g", prob.g, "eta_g_1", eta_11_vec[2],
		   "eta_g_2", eta_21_vec[2],"eta_g_3", eta_31_vec[2], 
                   marg.type, pair.type, visit.rate, sep="_")
  
  s1.out <- paste("s1.", as.character(file.tx), ".dat", sep="")
  if ( file.exists(s1.out) ) { unlink(s1.out) }
  
  s2.out <- paste("s2.", as.character(file.tx), ".dat", sep="")
  if ( file.exists(s2.out) ) { unlink(s2.out) }
  
  if(marg.type == "gene-common"){
    cat("nsim sec", 
        paste0("logalp1", 1:R1), "eta10",
        paste0("logalp2", 1:R2), "eta20",
        paste0("logalp3", 1:R3), "eta30", "eta_x", "eta_g",
        paste0("ASE.logalp1", 1:R1), "ASE.eta10",
        paste0("ASE.logalp2", 1:R2), "ASE.eta20",
        paste0("ASE.logalp3", 1:R3), "ASE.eta30", "ASE.eta_x", "ASE.eta_g",
        "\n", sep=" ", append=T, file=s1.out)
    
    cat("nsim sec", 
        "gam11", "rhodgr11", "gam22", "rhodgr22", "gam33", "rhodgr33", 
        "gam12", "rhodgr12", "gam13", "rhodgr13", "gam23", "rhodgr23", 
        paste0("ASE.", c("gam11", "rhodgr11", "gam22", "rhodgr22", "gam33", "rhodgr33", 
                         "gam12", "rhodgr12", "gam13", "rhodgr13", "gam23", "rhodgr23")),
        "\n", sep=" ", append=T, file=s2.out)
    
  }else if(marg.type=="type-specific"){
    
    cat("nsim sec", 
        paste0("logalp1", 1:R1), "eta10",
        paste0("logalp2", 1:R2), "eta20",
        paste0("logalp3", 1:R3), "eta30", "eta_x", paste0("eta_g_", 1:3),
        paste0("ASE.logalp1", 1:R1), "ASE.eta10",
        paste0("ASE.logalp2", 1:R2), "ASE.eta20",
        paste0("ASE.logalp3", 1:R3), "ASE.eta30", "ASE.eta_x", paste0("ASE.eta_g_", 1:3),
        "\n", sep=" ", append=T, file=s1.out)
    
    cat("nsim sec", 
        "gam11", "rhodgr11", "gam22", "rhodgr22", "gam33", "rhodgr33", 
        "gam12", "rhodgr12", "gam13", "rhodgr13", "gam23", "rhodgr23", 
        paste0("ASE.", c("gam11", "rhodgr11", "gam22", "rhodgr22", "gam33", "rhodgr33", 
                         "gam12", "rhodgr12", "gam13", "rhodgr13", "gam23", "rhodgr23")),
        "\n", sep=" ", append=T, file=s2.out)
  }
    
  
  for(m in 1:nsim){
    
    set.seed(m)
    
    t1 <- Sys.time()
    
    # --------- Data generation --------- #
    
    # generate XGmat
    XGmat = gen_XG.f(nsample, prob.x, prob.g)
    xgInd_11 = XGmat[,1]*XGmat[,2]; xgInd_10 = XGmat[,1]*(1-XGmat[,2])
    xgInd_01 = (1-XGmat[,1])*XGmat[,2]; xgInd_00 = (1-XGmat[,1])*(1-XGmat[,2])
    
    n_11 = sum(xgInd_11); n_10 = sum(xgInd_10); n_01 = sum(xgInd_01); n_00 = sum(xgInd_00)
    
    # generate Z_g_XG
    Z_11 = gen_Z_g_XG.f(n_11, pjoint_Z_g_11, reZmat)$genZmat
    Z_10 = gen_Z_g_XG.f(n_10, pjoint_Z_g_10, reZmat)$genZmat
    Z_01 = gen_Z_g_XG.f(n_01, pjoint_Z_g_01, reZmat)$genZmat
    Z_00 = gen_Z_g_XG.f(n_00, pjoint_Z_g_00, reZmat)$genZmat
    
    reXGmat = rbind(XGmat[xgInd_11==1,], XGmat[xgInd_10==1,],
                    XGmat[xgInd_01==1,], XGmat[xgInd_00==1,])
    trueZmat = rbind(Z_11, Z_10, Z_01, Z_00) 
    
    # generate T
    #system.time(trueTmat <- gen_T.f(K1, K2, K3, nsample, alp_1_vec, alp_2_vec, alp_3_vec, brks_1_vec, brks_2_vec, brks_3_vec,cor_mat))
    trueTmat <- gen_T_wei.f(K1, K2, K3, nsample, scale_1, kappa_1,  scale_2, kappa_2,
                            scale_3, kappa_3,cor_mat)
    
    # generate visitTime
    visitTime <- t(hpp.sim(visit.rate, visit.rate*2, num.sims=nsample))
    
    # generate indata_s1.list (each element corresponds to (j,k))
    indata_s1_1.list <- mclapply(1:K1, function(k){
      gen_indata_s1_jk.f(1, k, nsample, trueTmat[, k], trueZmat[,k], reXGmat, visitTime)
    }, mc.cores=2)
    
    indata_s1_2.list <- mclapply((1:K2), function(k){
      gen_indata_s1_jk.f(2, k, nsample, trueTmat[, K1+k], trueZmat[,K1+k], reXGmat, visitTime)
    }, mc.cores=2)
    
    indata_s1_3.list <- mclapply((1:K3), function(k){
      gen_indata_s1_jk.f(3, k, nsample, trueTmat[, K1+K2+k], trueZmat[, K1+K2+k], reXGmat, visitTime)
    }, mc.cores = 2)
    
    indata_s1.list <- c(indata_s1_1.list, indata_s1_2.list, indata_s1_3.list)
    
    # generate df_s2.list
    df_s2.list <- create_df_s2_jkjpkp.f(indata_s1.list, K1, K2, K3)
    
    cat("data generated!", "\n")
    
    # --------- two-stage estimation --------- #
    
    # stage I
    brks.list = list(brks_1_vec, brks_2_vec, brks_3_vec)
    
    t1.s1 <- Sys.time()
    res_s1 <- estParam_s1.f(indata_s1.list, brks.list, marg.type=marg.type)
    t2.s1 <- Sys.time()
    
    phi1_est <- res_s1$estimate
    
    cat("stage I done! spend ", round(difftime(t2.s1, t1.s1, units = "sec"),2), "\n")
    
    if(marg.type=="common"){
      
      phi1_est.list <- rep(list(phi1_est), 3)
      
    }else if(marg.type=="type-specific"){
      
      R = R1 + R2 + R3
      phi1_1 = phi1_est[c(1:(R1+1), R+3+(1:2))]
      phi1_2 = phi1_est[c(R1+1+(1:(R2+1)), R+3+c(1, 3))]
      phi1_3 = phi1_est[c(R1+R2+2+(1:(R3+1)), R+3+c(1, 4))]
      phi1_est.list <- list(phi1_1, phi1_2, phi1_3)
      
    }else if(marg.type=="gene-common"){
      
      #R1 = length(brks_1_vec)+1; R2 = length(brks_2_vec)+1; R3 = length(brks_3_vec)+1
      R = R1 + R2 + R3
      phi1_1 = phi1_est[c(1:(R1+1), R+3+(1:2))]
      phi1_2 = phi1_est[c(R1+1+(1:(R2+1)), R+3+c(1,2))]
      phi1_3 = phi1_est[c(R1+R2+2+(1:(R3+1)), R+3+c(1,2))]
      phi1_est.list <- list(phi1_1, phi1_2, phi1_3)
    }
    
    # stage II
    t1.s2 <- Sys.time()
    res_s2 <- estParam_s2.f(phi1_est.list, brks.list, df_s2.list, pair.type=pair.type)
    t2.s2 <- Sys.time()
      
    cat("stage II done! spend ", round(difftime(t2.s2, t1.s2, units = "sec"),2), "\n")
    
    # --------- inference --------- #
    # phi1_est.list
    # phi2_est.list
    phi2_est.list <- lapply(res_s2, function(x) x$estimate)
    Hess_phi2.list <- lapply(res_s2, function(x) x$hessian)
    
    indata_s1 = do.call("rbind", indata_s1.list)
    
    Hess_s1 <- res_s1$hessian
    #print(phi1_est)
    system.time(
    ts_ase <- ts_asvar.f(Hess_s1, Hess_phi2.list, 1e-05, nsample, 
                           phi1_est, brks.list, indata_s1, 
                           phi2_est.list, df_s2.list, marg.type, num.cores=6)
    )
    #est.ase.ts = sqrt(diag(ts_asvar))
    
    p1 = length(phi1_est)
    
    cat(m, round(difftime(t2.s1, t1.s1, units = "sec"),2), 
        res_s1$estimate, ts_ase$phi1,
        "\n", sep=" ", append=T, file=s1.out)
    
    cat(m, round(difftime(t2.s2, t1.s2, units = "sec"),2),
        as.vector(unlist(phi2_est.list)),ts_ase$phi2,
        "\n", sep=" ", append=T, file=s2.out)
    
    t2 <- Sys.time()
    
    print(t2-t1)
  }
} 


# test
J3_twoStage_script.f(nsim=100, nsample=2000, 
                     prob.g=0.05, prob.x=0.45,
                     K1=2, K2=6, K3=6, 
                     pi1=0.2, pi2=0.05, pi3=0.05,
                     eta_11_vec=c(0.2, 0), eta_21_vec=c(0.2, 0), eta_31_vec=c(0.2, 0),
                     SurvP.ctr1= 0.9, SurvP.ctr2=0.85, SurvP.ctr3=0.85,
                     kappa_1=1, kappa_2=1, kappa_3=1,
                     brks_1_vec=0.55, brks_2_vec=0.55, brks_3_vec=0.55,
                     or_11=1.2, or_22=1.2, or_33=1.2, or_12=1.05, or_13=1.05, or_23=1.05, 
                     tau_11=0.15, tau_22=0.15, tau_33=0.15, tau_12=0.05, tau_13=0.05, tau_23=0.05,
                     visit.rate=5,
                     marg.type = "gene-common", pair.type = "jjp")
