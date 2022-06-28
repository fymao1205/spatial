
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


nsample=2000; Jtype=2; KJjoint=2
gam2=c(0.5,0.5)
prob.g=0.1
pi.marg.j=c(0.65,0.75); visit.rate=20
pair.asso.type="jjp"; gene.effect.type="common"
delta=1e-03; 
s2.num.cores = ase.num.cores=2
  
  param.res.list <- param.gen.f(Jtype=Jtype, KJjoint=KJjoint, prob.g=prob.g,
                                pi.marg.j=pi.marg.j, gam2=gam2)
  
  FT.param.list <- param.res.list$FT.param.list
  Z.param.list <- param.res.list$Z.param.list
  X.param.list <- param.res.list$X.param.list
  G.param.list <- param.res.list$G.param.list
  
  Rj_vec <- unlist(lapply(1:length(FT.param.list$brks.list),function(j){length(FT.param.list$brks.list[[j]])+1}))
  
    m=1
    set.seed(m)
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
    
  t.end <- Sys.time()
  print(t.end)
  
