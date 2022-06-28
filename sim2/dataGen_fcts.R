# generate data :
# failure time T = (T_11, T_12, T_13, T_21, T_22, T_23)
# susceptibility Z = (Z_11, Z_12, Z_13, Z_21, Z_22, Z_23)
# baseline covariate X
# genetic marker G

expit.f <- function(x) {
  return( exp(x) / ( 1 + exp(x) ) )
}

logit.f <- function(x) {
  return( log( x / ( 1 - x) ) )
}

gen.data.f <- function(nsample, Jtype, KJjoint, 
                       FT.param.list, Z.param.list, 
                       X.param.list, G.param.list){
  
  # generate covariate X --------------------------------
  if(X.param.list$dist.name == "Bernoulli"){
    x.sim <- rbinom(nsample, 1, X.param.list$prob)
  }
  
  if(G.param.list$dist.name == "Bernoulli"){
    g.sim <- rbinom(nsample, 1, G.param.list$prob)
  }
  
  XG_11_index <- (x.sim==1)&(g.sim==1)
  XG_10_index <- (x.sim==1)&(g.sim==0)
  XG_01_index <- (x.sim==0)&(g.sim==1)
  XG_00_index <- (x.sim==0)&(g.sim==0)
  
  xg.ordered.sim <- rbind(cbind(x.sim,g.sim)[XG_11_index,],
                          cbind(x.sim,g.sim)[XG_10_index,],
                          cbind(x.sim,g.sim)[XG_01_index,],
                          cbind(x.sim,g.sim)[XG_00_index,])
  
  
  # generate susceptibility variable Zi -----------------
  
  p_marg <- Z.param.list$param.marg
  
  or_mat <- Z.param.list$or.mat
  
  # generate the p_mar of Zi given Xi and Gi
  
  p_mar_11 <- unlist(lapply(1:Jtype, function(j){
    xx_11 <- sum(p_marg[[j]])
    rep(expit.f(xx_11), KJjoint) # conditional P(Z_jk|X=1,G=1) instead of P(Z_jk, X=1, G=1)
  }))
  
  p_mar_10 <- unlist(lapply(1:Jtype, function(j){
    xx_10 <- sum(p_marg[[j]][-3])
    rep(expit.f(xx_10), KJjoint)
  }))
  
  p_mar_01 <- unlist(lapply(1:Jtype, function(j){
    
    xx_01 <- sum(p_marg[[j]][-2])
    rep(expit.f(xx_01), KJjoint)
    
  }))
  
  p_mar_00 <- unlist(lapply(1:Jtype, function(j){
    
    xx_00 <- sum(p_marg[[j]][-c(2,3)])
    rep(expit.f(xx_00), KJjoint)
    
  }))
  
  pre.p_mar_list <- list(p_mar_11, p_mar_10, p_mar_01, p_mar_00)
  p_joint.list <- mclapply(pre.p_mar_list, function(x){
    ObtainMultBinaryDist(odds = or_mat, marg.probs = x)
  }, mc.cores=2)
  
  p_joint_11 <- p_joint.list[[1]]
  p_joint_10 <- p_joint.list[[2]]
  p_joint_01 <- p_joint.list[[3]]
  p_joint_00 <- p_joint.list[[4]]
  
  nsample_11 <- sum(XG_11_index)
  nsample_10 <- sum(XG_10_index)
  nsample_01 <- sum(XG_01_index)
  nsample_00 <- nsample - nsample_11 - nsample_10 - nsample_01
  
  pre.z.sim.list <- list(list(ns=nsample_11, pj=p_joint_11),
                         list(ns=nsample_10, pj=p_joint_10),
                         list(ns=nsample_01, pj=p_joint_01),
                         list(ns=nsample_00, pj=p_joint_00))
  
  z.sim.list <- mclapply(pre.z.sim.list, function(x){
    RMultBinary(n = x$ns, mult.bin.dist = x$pj)$binary.sequences
  }, mc.cores = 4)
  
  z.sim <- do.call("rbind", z.sim.list)
  
  #*********************** using mclapply for faster computing ************************#
  #z.sim_11 <- RMultBinary(n = nsample_11, mult.bin.dist = p_joint_11)$binary.sequences#
  #z.sim_10 <- RMultBinary(n = nsample_10, mult.bin.dist = p_joint_10)$binary.sequences#
  #z.sim_01 <- RMultBinary(n = nsample_01, mult.bin.dist = p_joint_01)$binary.sequences#
  #z.sim_00 <- RMultBinary(n = nsample_00, mult.bin.dist = p_joint_00)$binary.sequences#
  #z.sim <- rbind(z.sim_11, z.sim_10, z.sim_01, z.sim_00)                              #
  #************************************************************************************#
  
  # generate failure time Ti given Zi=1 -----------------------------
  
  corr_mat <- FT.param.list$corr.mat
  mass_list <- FT.param.list$mass.list
  breaks_list <- FT.param.list$brks.list
  
  K_sum <- Jtype*KJjoint
  
  invPhi_ST <- mvrnorm(nsample, mu=rep(0, K_sum),Sigma=corr_mat)
  
  u.sim <- pnorm(invPhi_ST)
  
  ft.sim <- mclapply(1:Jtype, function(j){
    mass_j <- mass_list[[j]]
    breaks_j <- breaks_list[[j]]
    ft_j <- qpch(1-as.vector(u.sim[,(1:KJjoint)+KJjoint*(j-1)]), breaks_j, mass_j)
    ft_j.mat <- matrix(ft_j, nrow=nsample, ncol=KJjoint, byrow = 0)
    return(ft_j.mat)
  }, mc.cores=2)
  
  ft.sim.mat <- do.call("cbind", ft.sim)
  
  truez <- as.vector(t(z.sim))
  ft <- as.vector(t(ft.sim.mat))#as.vector(unlist(ft.sim))
  ft[truez==0] <- Inf
  
  return(list(XG=xg.ordered.sim, FT=ft, Z=truez))
  
}

gen.visit.f <- function(rate, njoints, ft){
  
  last.2visits.list <- lapply(1:njoints, function(i){
    g_vec <- hpp.sim(rate, rate*2)
    vt_vec <- c(g_vec[g_vec<=1], Inf)
    et <- ft[i]
    
    if(is.infinite(et)){
      l1 <- Inf 
      l2 <- vt_vec[max(which(vt_vec<et))]
    }else{
      l1 <- vt_vec[min(which(vt_vec>et))]
      l2 <- vt_vec[max(which(vt_vec<=et))]
    }
    
    l1 <- ifelse(et==0, vt_vec[2], l1)
    l2 <- ifelse(et==0, 0, l2)
    if(l2 >= l1){stop(print(c(l1,l2,et)))}
    
    c(l1,l2)
  })
  last.2visits <- do.call("rbind", last.2visits.list)
  return(last.2visits)
}

