
# ------------------------------------------------------------------------------------------
# script of the 1st set of simulation studies and results presented in Table III 
# ------------------------------------------------------------------------------------------

library(mipfp)
library(Rsolnp)
library(MASS)
library(NlcOptim)
library(eha) #
library(parallel)
library(mnormt)#
library(poisson)
source("dataGen_fcts.R")
source("TwoStage_est_fct.R")
source("prepare_param_sim.R")
source("TwoStage_inference_fct.R")


# script function for simulation studies, include:
# 1. data generation
# 2. estimation and inference
# 3. record results in .dat files
ts.script.f <- function(nsim, nsample, Jtype, KJjoint, gam2, prob.g,
                        pi.marg.j, visit.rate, pair.asso.type,gene.effect.type,
                        delta, s2.num.cores, ase.num.cores){
  
  param.res.list <- param.gen.f(Jtype=Jtype, KJjoint=KJjoint, prob.g=prob.g,
                                pi.marg.j=pi.marg.j, gam2=gam2)
  
  FT.param.list <- param.res.list$FT.param.list
  Z.param.list <- param.res.list$Z.param.list
  X.param.list <- param.res.list$X.param.list
  G.param.list <- param.res.list$G.param.list
  
  Rj_vec <- unlist(lapply(1:length(FT.param.list$brks.list),function(j){length(FT.param.list$brks.list[[j]])+1}))
  
  file.tx <- paste("2stg.n",nsample, "J",Jtype,"K", KJjoint, 
                   "gam21", gam2[1], "gam22", gam2[2],
                   "prob.g",prob.g, "pi1", pi.marg.j[1], "pi2", pi.marg.j[2],
                   gene.effect.type,
                   pair.asso.type, visit.rate,
                   sep="_")
  
  # file record for stage I estimation
  s1.out <- paste("s1.", as.character(file.tx), ".dat", sep="")
  if ( file.exists(s1.out) ) { unlink(s1.out) }
  
  # file record for stage II estimation
  s2.out <- paste("s2.", as.character(file.tx), ".dat", sep="")
  if ( file.exists(s2.out) ) { unlink(s2.out) }
  
  if(gene.effect.type=="common"){
    cat("nsim sec", paste0("logalp", rep(1:Jtype, Rj_vec), unlist(lapply(1:Jtype, function(j){
      1:Rj_vec[j] }))), paste0("gam0",1:Jtype), paste0("gam",1:Jtype),
      paste0("ASE.logalp", rep(1:Jtype, Rj_vec), unlist(lapply(1:Jtype, function(j){
        1:Rj_vec[j] }))), paste0("ASE.gam0",1:Jtype), paste0("ASE.gam",1:Jtype),
      "\n", sep=" ", append=T, file=s1.out)
  
    }else if(gene.effect.type=="block"){
    cat("nsim sec", paste0("logalp", rep(1:Jtype, Rj_vec), unlist(lapply(1:Jtype, function(j){
      1:Rj_vec[j] }))), paste0("gam0",1:Jtype), "gam1", paste0("gam2",1:Jtype),
      paste0("ASE.logalp", rep(1:Jtype, Rj_vec), unlist(lapply(1:Jtype, function(j){
        1:Rj_vec[j] }))), paste0("ASE.gam0",1:Jtype), "ASE.gam1", paste0("ASE.gam2",1:Jtype),
      "\n", sep=" ", append=T, file=s1.out)
  }  
  
  
  if(pair.asso.type=="jjp"){
    
    jjp.mat <- unique(t(combn(rep(1:Jtype,each=KJjoint),2)))
    cat("nsim sec", 
        paste0(rep(c("tan_pi/2sigma", "eta"),dim(jjp.mat)[1]), rep(paste0(jjp.mat[,1],jjp.mat[,2]),each=2)),
        paste0(rep(c("ASE.tan_pi/2sigma", "ASE.eta"),dim(jjp.mat)[1]), rep(paste0(jjp.mat[,1],jjp.mat[,2]),each=2)),
        "\n", sep=" ", append=T, file=s2.out)
  }else if(pair.asso.type=="common"){
    cat("nsim sec tan_pi/2sigma eta",
        "ASE.tan_pi/2sigma ASE.eta",
        "\n", sep=" ", append=T, file=s2.out)
  }else if(pair.asso.type=="full"){
    
    cat("nsim sec", paste0("tan_pi/2sigma", t(combn(1:(Jtype*KJjoint),2))[,1],t(combn(1:(Jtype*KJjoint),2))[,2]),
        paste0("eta", t(combn(1:(Jtype*KJjoint),2))[,1],t(combn(1:(Jtype*KJjoint),2))[,2]),
        "\n", sep=" ", append=T, file=s2.out)
  }
  
  
  for(sim in 1:nsim){
    
    set.seed(sim)
    t.begin <- Sys.time()
    t_pre_data_start <- Sys.time()
    # generate data 
    pre_data <- gen.data.f(nsample, Jtype, KJjoint, 
                           FT.param.list, Z.param.list, 
                           X.param.list, G.param.list)
    t_pre_data_stop <- Sys.time()
    
    print(t_pre_data_stop - t_pre_data_start)
    
    # prepare indata_s1
    indata_s1 <- get.data.s1.f(nsample, Jtype, KJjoint, XG=pre_data$XG, 
                               FT=pre_data$FT, Z=pre_data$Z, visit.rate=visit.rate)
    
    # stage I
    t_res_s1_start <- Sys.time()
    
    res_s1 <- est.param.s1.f(indata_s1, brks_list=FT.param.list$brks.list,
                             gene.effect.type=gene.effect.type)
    
    t_res_s1_stop <- Sys.time()
    print(t_res_s1_stop - t_res_s1_start)
    
    # prepare indata_s2
    system.time(indata_s2 <- get.data.s2.f(nsample, Jtype, KJjoint, outdata_s1=res_s1$outdata))
    
    # stage II
    t_res_s2_start <- Sys.time()
    res_s2 <- est.param.s2.f(Kj_vec=rep(KJjoint,Jtype), indata_s2=indata_s2, 
                             pair.asso.type=pair.asso.type, s2.num.cores=s2.num.cores)
    t_res_s2_stop <- Sys.time()
    print(t_res_s2_stop - t_res_s2_start)
    
    # ase calculation
    est.A11 <- res_s1$res$hessian
    est.A22 <- res_s2$res$hessian
    est.alpha.list <- list(res_s1$est.alp1_vec, res_s1$est.alp2_vec)
    if(gene.effect.type=="common"){
      est.gam <- rbind(c(res_s1$est.gam10, res_s1$est.gam1, res_s1$est.gam2[1]),
                       c(res_s1$est.gam20, res_s1$est.gam1, res_s1$est.gam2[1]))
    }else if(gene.effect.type=="block"){
      est.gam <- rbind(c(res_s1$est.gam10, res_s1$est.gam1, res_s1$est.gam2[1]),
                       c(res_s1$est.gam20, res_s1$est.gam1, res_s1$est.gam2[2]))
    }
   
    est.sig_vec <- res_s2$est.sigma_pw
    est.eta_vec <- res_s2$est.eta_pw

    est.asvar <- eval.asvar_ts.f(delta=1e-3, nsample=nsample, Kj_vec=rep(KJjoint,Jtype), est.A11=est.A11, est.A22=est.A22,
                                 indata_s1=indata_s1, est.alpha.list=est.alpha.list, brks_list=FT.param.list$brks.list, 
                                 est.gam=est.gam, gene.effect.type=gene.effect.type,
                                 indata_s2=indata_s2, est.sig_vec=est.sig_vec, est.eta_vec=est.eta_vec, 
                                 pair.asso.type=pair.asso.type, ase.num.cores=ase.num.cores)
    
    print(est.asvar)
    est.ase.ts <- sqrt(diag(est.asvar)) 
    
    cat(sim, round(difftime(t_res_s1_stop, t_res_s1_start, units = "sec"),2), 
        res_s1$res$estimate, est.ase.ts[(1:(sum(Rj_vec)+Jtype+2))],
        "\n", sep=" ", append=T, file=s1.out)
    
    cat(sim, round(difftime(t_res_s2_stop, t_res_s2_start, units = "sec"),2),
        res_s2$res$estimate,
        est.ase.ts[-(1:(sum(Rj_vec)+Jtype+2))],
        "\n", sep=" ", append=T, file=s2.out)
  }
  
  t.end <- Sys.time()
  print(t.end)
  
}


# test
ts.script.f(nsim=10, nsample=2000, Jtype=2, KJjoint=2, 
            gam2=c(0.5,0.5),
            prob.g=0.1,
            pi.marg.j=c(0.65,0.75), visit.rate=20,
            pair.asso.type="jjp", gene.effect.type="common",
            delta=1e-03, s2.num.cores = 2,ase.num.cores=2)

