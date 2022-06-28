
# Solve parameters for data generation

sol.gam_j0_vec.f <- function(pi_marg_j_vec, gam_1, gam_gj_vec, prob_x, prob_g){
  
  res.list <- sapply(1:length(pi_marg_j_vec), function(x){
    g <- function(gam_0){
      pi_marg_j_vec[x] - expit.f(gam_0+gam_1+gam_gj_vec[x])*prob_x*prob_g -
        expit.f(gam_0+gam_1)*prob_x*(1-prob_g) -
        expit.f(gam_0+gam_gj_vec[x])*(1-prob_x)*prob_g -
        expit.f(gam_0)*(1-prob_x)*(1-prob_g)
    }
    
    uniroot(g, lower=-10, upper=10)$root
    
  })
  
  res <- unlist(res.list)
  
  return(res)
  
}

solve.alp_j0 <- function(t.ctr_j, breaks_j, pre_alp_j, S.ctr_j){
  
  num <- -log(1-S.ctr_j)
  breaks.use <- (breaks_j <= t.ctr_j)
  brks.l <- c(0,breaks_j[breaks.use])
  brks.u <- c(brks.l[-1], t.ctr_j)
  len.brks <- brks.u -brks.l
  exp_pre_alp_j_serial <- exp(c(1:length(brks.l)-1)*pre_alp_j)
  denom <- sum(len.brks*exp_pre_alp_j_serial)
  mass_j0 <- num/denom
  
  return(mass_j0)
}

param.gen.f <- function(Jtype, KJjoint, prob.g,pi.marg.j, gam2){
  
  prob.x <- 0.5
  #prob.g <- prob.g # 0.25, 0.1
  
  # parameter list for X
  X.param.list <- list(dist.name="Bernoulli", prob=prob.x)
  
  # parameter list for G
  G.param.list <- list(dist.name="Bernoulli", prob=prob.g) #0.25, 0.1 
  
  # parameter list for Z
  # solve for gam_j0_vec = (gam_10, gam_20)
  pi_marg_j <- pi.marg.j # (0.1,0.25,0.5)
  
  OR11 <- 1.15; OR22 <- 1.1; OR12 <- 1.05 #1.05
  or.mat <- matrix(Inf, nrow=KJjoint*Jtype, ncol=KJjoint*Jtype)
  or.mat_11 <- matrix(OR11, nrow = KJjoint, ncol=KJjoint)
  or.mat_12 <- matrix(OR12, nrow = KJjoint, ncol=KJjoint)
  or.mat_22 <- matrix(OR22, nrow = KJjoint, ncol=KJjoint)
  or.mat <- rbind(cbind(or.mat_11, or.mat_12), cbind(or.mat_12, or.mat_22))
  diag(or.mat) <- Inf
  gam_1 <- 0.45
  gam_gj_vec <- gam2
  
  gam_j0_vec <- sol.gam_j0_vec.f(pi_marg_j_vec=pi_marg_j, gam_1=gam_1, gam_gj_vec=gam_gj_vec, 
                                 prob_x=prob.x, prob_g=prob.g)
  
  Z.param.list <- list(param.marg=list(c(gam_j0_vec[1], gam_1, gam_gj_vec[1]),
                                       c(gam_j0_vec[2], gam_1, gam_gj_vec[2])),
                       or.mat=or.mat)
  
  # parameter list for FT
  
  # kendall's tau
  tau1 <- 0.1; tau2 <- 0.15; tau12 <- 0 #asin(2/pi)
  corr_11 <- pi/2*sin(tau1)
  corr_22 <- pi/2*sin(tau2)
  corr_12 <- pi/2*sin(tau12)
  Bv_corr_U_11 <- matrix(corr_11,nrow = KJjoint, ncol=KJjoint)
  diag(Bv_corr_U_11) <- 1
  Bv_corr_U_22 <- matrix(corr_22,nrow = KJjoint, ncol=KJjoint)
  diag(Bv_corr_U_22) <- 1
  Bv_corr_U_12 <- matrix(corr_12,nrow = KJjoint, ncol=KJjoint)
  Bv_corr_U_r1 <- cbind(Bv_corr_U_11, Bv_corr_U_12)
  Bv_corr_U_r2 <- cbind(Bv_corr_U_12, Bv_corr_U_22)
  Bv_corr_U_mat <- rbind(Bv_corr_U_r1, Bv_corr_U_r2)
  
  pre_alp_1 <- log(1)
  pre_alp_2 <- log(1.1)
  Rj <- 2; breaks.vec <- c(0.55) # (0,0.35), (0.35,0.70), (0.70, Inf)
  t.ctr <- 1
  S1.ctr <- 0.85
  S2.ctr <- 0.9
  
  # using constraints to solve for 'intercept'
  alp_10 <- solve.alp_j0(1,breaks.vec, pre_alp_1, S1.ctr) # solve for \rho_1
  alp_20 <- solve.alp_j0(1,breaks.vec, pre_alp_2, S2.ctr) # solve for \rho_2
  
  # pwc mass for joint type j, j=1,2,3
  alp_1_vec <- alp_10*exp(pre_alp_1*c(1:(Rj)-1))
  alp_2_vec <- alp_20*exp(pre_alp_2*c(1:(Rj)-1))
  
  FT.param.list <- list(corr.mat=Bv_corr_U_mat,
                        mass.list=list(alp_1_vec,alp_2_vec),
                        brks.list=list(breaks.vec, breaks.vec))
  
  return(list(X.param.list=X.param.list,
         G.param.list=G.param.list,
         Z.param.list=Z.param.list,
         FT.param.list=FT.param.list))

}



