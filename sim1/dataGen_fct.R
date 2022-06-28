
# 0. J=3: K1, K2, K3; K= K1+K2+K3, total # of joints
# 1. prepare parameters
# 2. generate data: X, Z, T, l, r
# 3. prepare indata_s1 for stage I estimation
# 4. prepare indata_s2 for stage II estimation (may be using Rcpp)

expit.f <- function(x) {
  return( exp(x) / ( 1 + exp(x) ) )
}

logit.f <- function(x) {
  return( log( x / ( 1 - x) ) )
}


# --------- 1. solve for necessary params
# --- 1.1. marginal models of Z_jk: P(Z_jk=1|X) = expit(eta0 + eta_1j*X); 
# -------  given eta_1j <2 x 1 vector>, solve eta0 <scalar> by setting E_X[P(Z_jk=1|X)] = p_j

sol.eta_0.f <- function(p_j, eta_1j, prob_x){
  
  g <- function(eta_0){
    p_j - expit.f(eta_0+sum(eta_1j))*prob_x[1]*prob_x[2] -
      expit.f(eta_0+eta_1j[1])*prob_x[1]*(1-prob_x[2]) -
      expit.f(eta_0+eta_1j[2])*prob_x[2]*(1-prob_x[1]) -
      expit.f(eta_0)*(1-prob_x[1])*(1-prob_x[2])
  }
  
  res = uniroot(g, lower=-10, upper=10)$root
  
  return(res)
}

# --- 1.2. marginal models for T_jk|Z_jk=1: Surv_j.ctr is actually CDF_j.ctr
# P(T_jk<=t|Z_jk=1) = pwc hazard with cuts brks_j and levels alp_j
# -------  alp_j1 = alp_0*exp(kappa_j); set kappa_j and solve alp_0
sol.alp_0.f <- function(Surv_j.ctr, t0_j.ctr, brks_j, kappa_j){
  
  g <- function(alp_0){
    
    levels = alp_0*exp(kappa_j*((1:(length(brks_j)+1))-1))
    
    survP = ppch(t0_j.ctr, cuts=brks_j, levels=levels)
    
    survP - Surv_j.ctr
  }
  
  res = uniroot(g, lower=-10, upper=10)$root
  
  return(res)
}

# ------- this is what used for data generation (as qpch() is too slow ...) 
# ------- for weibull dist. (scale, shape); in pweibull(scale=1/scale, shape=shape)
sol.scale.f <- function(Surv_j.ctr, t0_j.ctr, shape){
  
  g <- function(scale){
    
    survP = pweibull(t0_j.ctr, scale=1/scale, shape=shape)
    
    survP - Surv_j.ctr
  }
  
  res = uniroot(g, lower=0, upper=10)$root
  
  return(res)
  
}

# test; result: the same ...
#sol.alp_0.f(0.85, 1, 0.55, log(1))
#solve.alp_j0(1, 0.55, log(1), 0.85)

# -------- 2. generate X, Z|X, T|Z ....
# --- 2.1. generate covarate X and gene G, binary and independent 
gen_XG.f <- function(n, prob_x, prob_g){
  
  x = rbinom(n, 1, prob_x)
  
  g = rbinom(n, 1, prob_g)
  
  res = cbind(x, g)
  
  return(res)
}

# --- 2.2. generate Z|(1,1), Z|(1,0), Z|(0,1), Z|(0,0)
# ------- 2.2.1. calculate the joint dist. for Z|(x,g)
J3_sol_Pm.f <- function(Binary_vec, pi_marg_vec){
  
  m = sum(Binary_vec)
  
  if(m<3){ stop("sum(Z) needs to be more than 3!")}
  res = prod(pi_marg_vec[Binary_vec==1])*prod(1-pi_marg_vec[Binary_vec==0])
  #res = pi1^m1*(1-pi1)^(2-m1)*pi2^m2*(1-pi2)^(2-m2)*pi3^m3*(1-pi3)^(2-m3)
  return(res)
}

J3_sol_P2.f <- function(Binary_vec, pi_marg_vec, pi_pair_mat){
  
  m = sum(Binary_vec)
  
  if(m!=2){ stop("sum(Z) needs to be 2!")}
  
  res2 = prod(pi_marg_vec[Binary_vec==1])*prod(1-pi_marg_vec[Binary_vec==0])
  
  pairInd = which(Binary_vec==1)
  res = pi_pair_mat[pairInd[1], pairInd[2]] - prod(pi_marg_vec[Binary_vec==1]) + res2
  
  return(res)
}

J3_sol_P1_1.f <- function(K1, K2, K3, pi11, pi12, pi13, pi1, pi2, pi3){
  
  res = (pi1*(1-pi1)^(K1-1)*(1-pi2)^(K2)*(1-pi3)^(K3) - (K1-1)*pi11 - K2*pi12 - K3*pi13 
         + pi1*( (K1-1)*pi1 + K2*pi2 +K3*pi3 ))
  
  return(res)
}

J3_sol_P1_2.f <- function(K1, K2, K3, pi12, pi22, pi23,pi1, pi2, pi3){
  
  res = (pi2*(1-pi2)^(K2-1)*(1-pi1)^(K1)*(1-pi3)^(K3) - (K2-1)*pi22 - K1*pi12 - K3*pi23 
         + pi2*( (K2-1)*pi2 + K1*pi1 +K3*pi3 ))
  
  return(res)
}

J3_sol_P1_3.f <- function(K1, K2, K3, pi13, pi23, pi33, pi1, pi2, pi3){
  
  #res <- pi3*(1-pi3)*((1-pi1)*(1-pi2))^2  - pi33 - pi13*2 - pi23*2 + pi3*(pi1*2+pi2*2+pi3)
  
  res = (pi3*(1-pi3)^(K3-1)*(1-pi1)^(K1)*(1-pi2)^(K2) - (K3-1)*pi33 - K1*pi13 - K2*pi23 
         + pi3*( (K3-1)*pi3 + K1*pi1 +K2*pi2 ))
  
  return(res)
}

create_pi_pair_mat.f <- function(K1, K2, K3, pi1, pi2, pi3, or11, or22, or33, or12, or13, or23){
  
  pi11 = pi11_jkjpkp(pi1, pi1, log(or11)); pi12 =  pi11_jkjpkp(pi1, pi2, log(or12)); pi13 =  pi11_jkjpkp(pi1, pi3, log(or13))
  pi22 = pi11_jkjpkp(pi2, pi2, log(or22)); pi23 =  pi11_jkjpkp(pi2, pi3, log(or23)); pi33 =  pi11_jkjpkp(pi3, pi3, log(or33))
  
  K = K1 + K2 + K3
  pi_pair_mat <- matrix(NA, K, K)
  or_mat <- matrix(or33, K, K)
  or_mat[1:K1, 1:K1] <- or11
  or_mat[(K1+1):(K1+K2), (K1+1):(K1+K2)] <- or22
  or_mat[1:K1, (K1+1):(K1+K2)] <- or12
  or_mat[1:K1, (K1+K2+1):(K1+K2+K3)] <- or13
  or_mat[(K1+1):(K1+K2), (K1+K2+1):(K1+K2+K3)] <- or23
  diag(or_mat) = NA
  
  pi_pair_mat <- matrix(pi33, K, K)
  pi_pair_mat[1:K1, 1:K1] <- pi11
  pi_pair_mat[(K1+1):(K1+K2), (K1+1):(K1+K2)] <- pi22
  pi_pair_mat[1:K1, (K1+1):(K1+K2)] <- pi12
  pi_pair_mat[1:K1, (K1+K2+1):(K1+K2+K3)] <- pi13
  pi_pair_mat[(K1+1):(K1+K2), (K1+K2+1):(K1+K2+K3)] <- pi23
  diag(pi_pair_mat) = NA
  
  for(i in 1:(K-1)){
    
    for(j in (i+1):K){
      
      or_mat[j, i] = or_mat[i,j] 
      pi_pair_mat[j, i] = pi_pair_mat[i,j] 
      
    }
  }
  
  
  res <- list(pi_pair_mat=pi_pair_mat,
              or_mat = or_mat,
              pi11 = pi11,
              pi12 = pi12,
              pi13 = pi13,
              pi22 = pi22,
              pi23 = pi23,
              pi33 = pi33)
  
  return(res)
}

J3_sol_pi_marg.f <- function(Zmat, K1, K2, K3, pi1, pi2, pi3, pi_pair_mat,
                             pi11, pi22, pi33, pi12, pi13, pi23){
  
  pi_marg_vec = c(rep(pi1,K1), rep(pi2, K2), rep(pi3, K3))
  
  #P6 = (pi1*pi2*pi3)^2
  #Z6mat_5 = Z6mat[rowSums(Z6mat)==5,]
  #P5 = apply(Z6mat_5, 1, K6J3_sol_Pm.f, pi_marg_vec=pi_marg_vec)
  #Z6mat_4 = Z6mat[rowSums(Z6mat)==4,]
  #P4 = apply(Z6mat_4, 1, K6J3_sol_Pm.f, pi_marg_vec=pi_marg_vec)
  Zmat_3 = Zmat[rowSums(Zmat)>=3,]
  P3 = apply(Zmat_3, 1, J3_sol_Pm.f, pi_marg_vec=pi_marg_vec)
  Zmat_2 = Zmat[rowSums(Zmat)==2,]
  P2 = apply(Zmat_2, 1, J3_sol_P2.f, pi_marg_vec=pi_marg_vec, pi_pair_mat = pi_pair_mat)
  
  #pi_pair_vec <- unique(pi_pair_mat[upper.tri(pi_pair_mat)])
  #pi11 <- pi_pair_vec[1]; 
  #pi12 <- pi_pair_vec[2]; pi22 <- pi_pair_vec[3]; 
  #pi23 <- pi_pair_vec[4]; pi13 <- pi_pair_vec[5]; pi33 <- pi_pair_vec[6]
  
  P1_1 = J3_sol_P1_1.f(K1, K2, K3, pi11, pi12, pi13, pi1, pi2, pi3)
  P1_2 = J3_sol_P1_2.f(K1, K2, K3, pi12, pi22, pi23, pi1, pi2, pi3)
  P1_3 = J3_sol_P1_3.f(K1, K2, K3, pi13, pi23, pi33, pi1, pi2, pi3)
  
  P0 = 1- sum(P3) - sum(P2) - (K1*P1_1+ K2*P1_2+ K3*P1_3)
  
  res <- c(P0, rep(P1_1,K1), rep(P1_2,K2), rep(P1_3, K3), P2, P3)
  
  return(res)
}


J3_Pjoint_Z_g_XG.f <- function(K1, K2, K3, eta_1_vec, eta_2_vec, eta_3_vec, x, g, 
                               pi_pair_mat,  pi11, pi22, pi33, pi12, pi13, pi23){
                               #or11, or22, or33, or12, or13, or23
                              
  # calculate marg prob.s pi1(x,g), pi2(x,g), pi3(x,g)
  pi1 = expit.f(sum(eta_1_vec*c(1,x,g))) 
  pi2 = expit.f(sum(eta_2_vec*c(1,x,g))) 
  pi3 = expit.f(sum(eta_3_vec*c(1,x,g))) 
  
  K= K1+K2+K3
  l <- rep(list(0:1), K)
  Zmat <- expand.grid(l)
  pjoint_vec = J3_sol_pi_marg.f(Zmat, K1, K2, K3, pi1, pi2, pi3, pi_pair_mat, 
                                pi11, pi22, pi33, pi12, pi13, pi23)
  
  return(pjoint_vec)
}

# ------- 2.2.2. generate Z|(x,g) given P(Z|x,g)
gen_Z_g_XG.f <- function(n_xg, pjoint_vec, reZmat){
  
  # generate random numbers from unif(0,1)
  u_vec = runif(n_xg)
  
  cutp_vec = cumsum(pjoint_vec)
  
  categ_vec = sapply(u_vec, function(x){
     sum(x>=cutp_vec)+1
  })
  
  genZmat = reZmat[categ_vec,]
  
  res = list(genZmat = genZmat, 
             categ_vec = categ_vec)
  
  return(res)
}

# --- 2.3. generate T
create_cor_mat.f <- function(K1, K2, K3, tau11, tau22, tau33, tau12, tau13, tau23){
  
  cor11 <- pi/2*sin(tau11); cor22 <- pi/2*sin(tau22); cor33 <- pi/2*sin(tau33)
  cor12 <- pi/2*sin(tau12); cor13 <- pi/2*sin(tau13); cor23 <- pi/2*sin(tau23)
  
  K = K1 + K2 + K3
  cor_mat <- matrix(NA, K, K)
  cor_mat <- matrix(cor33, K, K)
  cor_mat[1:K1, 1:K1] <- cor11
  cor_mat[(K1+1):(K1+K2), (K1+1):(K1+K2)] <- cor22
  cor_mat[1:K1, (K1+1):(K1+K2)] <- cor12
  cor_mat[1:K1, (K1+K2+1):K] <- cor13
  cor_mat[(K1+1):(K1+K2), (K1+K2+1):K] <- cor23
  diag(cor_mat) = 1
  
  for(i in 1:(K-1)){
    
    for(j in (i+1):K){
      
      cor_mat[j, i] = cor_mat[i,j] 
    }
  }
  
  return(cor_mat)
}

# ------- 2.3.1. generate true T
# slow .......
# Rcpp is still slow ... b/c the established R function qpch() is slow ...
gen_T.f <- function(K1, K2, K3, n, #Z1, Z2, Z3, 
                    alp_1_vec, alp_2_vec, alp_3_vec, brks_1, brks_2, brks_3,
                    cor_mat){
  
  K = K1+K2+K3
  
  invPhi_ST <- mvrnorm(n, mu=rep(0, K),Sigma=cor_mat)
  
  u.sim <- pnorm(invPhi_ST)
  
  # this is slow ...; no covariate is involved
  #  Rcpp is still slow ..., try mclapply, a little bit improved
  #microbenchmark(mclapply(1-as.vector(u.sim[, 1:K1]), qpch, cuts=brks_1, levels=alp_1_vec, mc.cores = 3), 
  #               qpch(1-as.vector(u.sim[, 1:K1]), brks_1, alp_1_vec), times=2)
  ft_1 = mclapply(1-as.vector(u.sim[, 1:K1]), qpch, cuts=brks_1, levels=alp_1_vec, mc.cores = 3)
  ft_2 = mclapply(1-as.vector(u.sim[, (K1+1):(K1+K2)]), qpch, cuts=brks_2, levels=alp_2_vec, mc.cores = 3)
  ft_3 = mclapply(1-as.vector(u.sim[, (K1+K2+1):K]), qpch, cuts=brks_3, levels=alp_3_vec, mc.cores = 3)
  
  #ft_1 = qpch(1-as.vector(u.sim[, 1:K1]), brks_1, alp_1_vec)
  #ft_2 = qpch(1-as.vector(u.sim[, (K1+1):(K1+K2)]), brks_2, alp_2_vec)
  #ft_3 = qpch(1-as.vector(u.sim[, (K1+K2+1):K]), brks_3, alp_3_vec)
  
  #truez_1 <- as.vector(Z1)
  ft_1_vec <- as.vector(unlist(ft_1))#as.vector(unlist(ft.sim))
  #ft_1_vec[truez_1==0] <- Inf
  
  #truez_2 <- as.vector(Z2)
  ft_2_vec <- as.vector(unlist(ft_2))#as.vector(unlist(ft.sim))
  #ft_2_vec[truez_2==0] <- Inf
  
  #truez_3 <- as.vector(Z3)
  ft_3_vec <- as.vector(unlist(ft_3))#as.vector(unlist(ft.sim))
  #ft_3_vec[truez_3==0] <- Inf
  
  #ft_res <- list(ft_1=ft_1, ft_2=ft_2, ft_3=ft_3)
  ft_1_mat <- matrix(ft_1_vec, ncol=K1, byrow=0) 
  ft_2_mat <- matrix(ft_2_vec, ncol=K2, byrow=0) 
  ft_3_mat <- matrix(ft_3_vec, ncol=K3, byrow=0) 
  
  ft_res <- cbind(ft_1_mat, ft_2_mat, ft_3_mat)
  
  return(ft_res)
}

gen_T_wei.f <- function(K1, K2, K3, n, scale_1, shape_1,  scale_2, shape_2, 
                         scale_3, shape_3, cor_mat){
  
  K = K1+K2+K3
  
  invPhi_ST <- mvrnorm(n, mu=rep(0, K), Sigma=cor_mat)
  
  u.sim <- pnorm(invPhi_ST)
  #ft_1 = mclapply(1-as.vector(u.sim[, 1:K1]), qpch, cuts=brks_1, levels=alp_1_vec, mc.cores = 3)
  ft_1 = qweibull(as.vector(u.sim[, 1:K1]), shape=shape_1, scale=1/scale_1, lower.tail = FALSE)
  ft_2 = qweibull(as.vector(u.sim[, K1+(1:K2)]), shape=shape_2, scale=1/scale_2, lower.tail = FALSE)
  ft_3 = qweibull(as.vector(u.sim[, K1+K2+(1:K3)]), shape=shape_3, scale=1/scale_3, lower.tail = FALSE)
  
  ft_1_mat <- matrix(ft_1, ncol=K1, byrow=0) 
  ft_2_mat <- matrix(ft_2, ncol=K2, byrow=0) 
  ft_3_mat <- matrix(ft_3, ncol=K3, byrow=0) 
  
  ft_res <- cbind(ft_1_mat, ft_2_mat, ft_3_mat)
  
  return(ft_res)
  
}

# --- 2.4. generate visit time and then (l, r); individual-level visit time
# -------- generate last two visit time for each joint
# slow .... think about Rcpp it ....
# pre operate: vt_vec <- c(visit_vec[visit_vec<=1], Inf) 
# Rcpp version <gen_last2visit()> much faster !!!
gen.last2visit.f <- function(visit_mat, n, ft){
  
  last.2visits.list <- mclapply(1:n, function(i){
    g_vec <- visit_mat[i,] #hpp.sim(rate, rate*2)
    vt_vec <- c(g_vec[g_vec<=1], Inf)
    et <- ft[i]
    
    if(is.infinite(et)){
      l1 <- Inf #vt_vec[length(vt_vec)-1]
      l2 <- vt_vec[max(which(vt_vec<et))]
    }else{
      l1 <- vt_vec[min(which(vt_vec>et))]
      l2 <- vt_vec[max(which(vt_vec<=et))]
    }
    
    #c(l2,l1)
    #interval.l <- ifelse(ft==0, 0, interval.l)
    #interval.u <- ifelse(ft==0, 0.1, interval.u)
    l1 <- ifelse(et==0, vt_vec[2], l1)
    l2 <- ifelse(et==0, 0, l2)
    if(l2 >= l1){stop(print(c(l1,l2,et)))}
    
    # c(upperbound, lowerbound)
    c(l1,l2)
  }, mc.cores=2)
  
  last.2visits <- do.call("rbind", last.2visits.list)
  return(last.2visits)
  
}

# --------- 3. prepare indata_s1 for stage I estimation; 
# --- 3.1. generate dataframe for each joint (j,k), nrow = n
# visitTime <- hpp.sim(rate, rate*2, num.sims=n); common for each (j,k)
# reXGmat: reorganized XGmat: (1,1), (1,0), (0,1), (0,0); correspond to trueZvec_jk
gen_indata_s1_jk.f <- function(j, k, n, trueTvec_jk, trueZvec_jk,
                               #trueZvec_11_jk, trueZvec_10_jk,
                               #trueZvec_01_jk, trueZvec_00_jk,
                               reXGmat, visitTime){
  
  joint = rep(k, n); type = rep(j, n); id = 1:n
  
  obsTvec_jk <- trueTvec_jk
  obsTvec_jk <- ifelse(trueZvec_jk==1, trueTvec_jk, Inf)
  last2visit <- gen_last2visit(visitTime, n, obsTvec_jk)
  interval.u <- last2visit[,1]
  interval.l <- last2visit[,2]
  
  indata <- as.data.frame(cbind(joint, type, id, reXGmat, trueTvec_jk, trueZvec_jk, interval.l, interval.u))
  colnames(indata) <- c("joint", "type", "id",  "x", "g", "truet", "truez", "l", "r")
  return(indata)
    
}

# --- 3.2. generate a list of dataframes (corresponding to each joint (j,k))
# -------- in script: 3.2.1. prepare params for marginal models
# --------            3.2.2. prepare pairwise params
# --------            3.2.3. generate X, G, T, Z, visitTime
# --------            3.2.4. create a list of dataframes using gen_indata_s1_jk.f()
# --------                   1:K1, K1+1:K1+K2, K1+K2:K

# --------- 4. prepare indata_s2 for stage II estimation;
# --- 4.1. generate dataframe for each pair of joints (j,k) and (j', k')
create_df_s2_jkjpkp.f <- function(indata_s1.list, K1, K2, K3){
  
  K = K1+K2+K3
  
  kkp_11 = t(combn(1:K1, 2))
  kkp_22 = t(combn( K1 + (1:K2), 2))
  kkp_33 = t(combn( K1 + K2+ (1:K3), 2))
  
  kkp_12 = expand.grid(1:K1, (K1+1):(K1+K2))
  kkp_13 = expand.grid(1:K1, (K1+K2+1):K)
  kkp_23 = expand.grid((K1+1):(K1+K2), (K1+K2+1):K)
  
  # kkp_mat <- rbind(kkp_11, kkp_22, kkp_33, kkp_12, kkp_13, kkp_23)
  
  df_s2_11.list = apply(kkp_11, 1, function(x){
     res <- as.data.frame(cbind(indata_s1.list[[ x[1] ]], 
                                indata_s1.list[[ x[2] ]][, c("joint", "type", "truet", "truez", "l", "r")]))
     colnames(res) <- c("joint_k", "type_j", "id", "x", "g", "truet_j", "truez_j", "l_j", "r_j",
                        "joint_kp", "type_jp", "truet_jp", "truez_jp", "l_jp", "r_jp")
     res
  })
  
  df_s2_22.list = apply(kkp_22, 1, function(x){
    res <- as.data.frame(cbind(indata_s1.list[[ x[1] ]], 
                               indata_s1.list[[ x[2] ]][, c("joint", "type", "truet", "truez", "l", "r")]))
    colnames(res) <- c("joint_k", "type_j", "id", "x", "g", "truet_j", "truez_j", "l_j", "r_j",
                       "joint_kp", "type_jp", "truet_jp", "truez_jp", "l_jp", "r_jp")
    res
  })
  
  df_s2_33.list = apply(kkp_33, 1, function(x){
    res <- as.data.frame(cbind(indata_s1.list[[ x[1] ]], 
                               indata_s1.list[[ x[2] ]][, c("joint", "type", "truet", "truez", "l", "r")]))
    colnames(res) <- c("joint_k", "type_j", "id", "x", "g", "truet_j", "truez_j", "l_j", "r_j",
                       "joint_kp", "type_jp", "truet_jp", "truez_jp", "l_jp", "r_jp")
    res
  })
  
  df_s2_12.list = apply(kkp_12, 1, function(x){
    res <- as.data.frame(cbind(indata_s1.list[[ x[1] ]], 
                               indata_s1.list[[ x[2] ]][, c("joint", "type", "truet", "truez", "l", "r")]))
    colnames(res) <- c("joint_k", "type_j", "id", "x", "g", "truet_j", "truez_j", "l_j", "r_j",
                       "joint_kp", "type_jp", "truet_jp", "truez_jp", "l_jp", "r_jp")
    res
  })
  
  df_s2_13.list = apply(kkp_13, 1, function(x){
    res <- as.data.frame(cbind(indata_s1.list[[ x[1] ]], 
                               indata_s1.list[[ x[2] ]][, c("joint", "type", "truet", "truez", "l", "r")]))
    colnames(res) <- c("joint_k", "type_j", "id", "x", "g", "truet_j", "truez_j", "l_j", "r_j",
                       "joint_kp", "type_jp", "truet_jp", "truez_jp", "l_jp", "r_jp")
    res
  })
  
  df_s2_23.list = apply(kkp_23, 1, function(x){
    res <- as.data.frame(cbind(indata_s1.list[[ x[1] ]], 
                               indata_s1.list[[ x[2] ]][, c("joint", "type", "truet", "truez", "l", "r")]))
    colnames(res) <- c("joint_k", "type_j", "id", "x", "g", "truet_j", "truez_j", "l_j", "r_j",
                       "joint_kp", "type_jp", "truet_jp", "truez_jp", "l_jp", "r_jp")
    res
  })
  
  df_s2_11 = do.call("rbind", df_s2_11.list)
  df_s2_22 = do.call("rbind", df_s2_22.list)
  df_s2_33 = do.call("rbind", df_s2_33.list)
  df_s2_12 = do.call("rbind", df_s2_12.list)
  df_s2_13 = do.call("rbind", df_s2_13.list)
  df_s2_23 = do.call("rbind", df_s2_23.list)
  
  df_s2.list = list(df_s2_11, df_s2_22, df_s2_33,
                    df_s2_12, df_s2_13, df_s2_23)
  
  return(df_s2.list)
}

# (it seems that using Rcpp is not too slow to prepare indata_s2_jkjpkp within log lik evaluation...)
# (so hold on ...)
# --- 4.2. prepare dataframe for insdata_s2_jjp 
# -------- using Rcpp function to generate indata_s2_jkjpkp


#microbenchmark(gen.last2visit.f(visitTime, n, trueTmat[,1]),
#               gen_last2visit(visitTime, n, trueTmat[,1]), times=2)


# --------- 5. check for generated data
# ---   prob.g, prob.x, 
# ---   pi1, pi2, pi3, or11, or22, or33, or12, or13, or23
# ---   SurvP.ctr1, SurvP.ctr2, SurvP.ctr3
dataGen_check <- function(nsim, nsample, 
                          prob.g, prob.x,
                          K1, K2, K3, 
                          pi1, pi2, pi3,
                          eta_11_vec, eta_21_vec, eta_31_vec,
                          SurvP.ctr1, SurvP.ctr2, SurvP.ctr3,
                          kappa_1, kappa_2, kappa_3,
                          or_11, or_22, or_33, or_12, or_13, or_23, 
                          tau_11, tau_22, tau_33, tau_12, tau_13, tau_23,
                          visit.rate){
  
  # prepare stage I params: alp_j_vec, eta_j_vec
  eta_10 = sol.eta_0.f(pi1, eta_11_vec, c(prob.x, prob.g))
  eta_20 = sol.eta_0.f(pi2, eta_21_vec, c(prob.x, prob.g))
  eta_30 = sol.eta_0.f(pi3, eta_31_vec, c(prob.x, prob.g))
  
  #alp_10 = sol.alp_0.f(SurvP.ctr1, 1, brks_1_vec, kappa_1)
  #alp_20 = sol.alp_0.f(SurvP.ctr2, 1, brks_2_vec, kappa_2)
  #alp_30 = sol.alp_0.f(SurvP.ctr3, 1, brks_3_vec, kappa_3)
  scale_1 = sol.scale.f(SurvP.ctr1, 1, kappa_1)
  scale_2 = sol.scale.f(SurvP.ctr2, 1, kappa_2)
  scale_3 = sol.scale.f(SurvP.ctr3, 1, kappa_3)
  
  eta_1_vec = c(eta_10, eta_11_vec); eta_2_vec = c(eta_20, eta_21_vec); eta_3_vec = c(eta_30, eta_31_vec)
  #alp_1_vec = alp_10*exp(kappa_1*(0:length(brks_1_vec)))
  #alp_2_vec = alp_20*exp(kappa_2*(0:length(brks_2_vec)))
  #alp_3_vec = alp_10*exp(kappa_3*(0:length(brks_3_vec)))
  
  # prepare stage II params: pi_pair_mat, cor_mat, Pjoint_Z_g_XG
  pi_pair_mat.res = create_pi_pair_mat.f(K1, K2, K3, or_11, or_22, or_33, or_12, or_13, or_23)
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
  if(any(pjoint_Z_g_01<0)){stop("pjoint_Z_g_11 is not legistimate !")}
  
  pjoint_Z_g_00 = J3_Pjoint_Z_g_XG.f(K1, K2, K3, eta_1_vec, eta_2_vec, eta_3_vec, 0, 0, pi_pair_mat,
                                     pi_pair_mat.res$pi11, pi_pair_mat.res$pi22, 
                                     pi_pair_mat.res$pi33, pi_pair_mat.res$pi12, 
                                     pi_pair_mat.res$pi13, pi_pair_mat.res$pi23)
  if(any(pjoint_Z_g_01<0)){stop("pjoint_Z_g_11 is not legistimate !")}
  
  res.mat <- matrix(NA, nsim, 28)
  for(i in 1:nsim){
    
    # generate XGmat
    XGmat = gen_XG.f(nsample, prob.x, prob.g)
    xgInd_11 = XGmat[,1]*XGmat[,2]; xgInd_10 = XGmat[,1]*(1-XGmat[,2])
    xgInd_01 = (1-XGmat[,1])*XGmat[,2]; xgInd_00 = (1-XGmat[,1])*(1-XGmat[,2])
    
    n_11 = sum(xgInd_11); n_10 = sum(xgInd_10); n_01 = sum(xgInd_01); n_00 = sum(xgInd_00)
    
    K= K1+K2+K3
    l <- rep(list(0:1), K)
    Zmat <- expand.grid(l)
    reZmat = rbind(Zmat[rowSums(Zmat)<=1,],
                   Zmat[rowSums(Zmat)==2,],
                   Zmat[rowSums(Zmat)>=3,])
    
    # generate Z_g_XG
    Z_11 = gen_Z_g_XG.f(n_11, pjoint_Z_g_11, reZmat)$genZmat
    Z_10 = gen_Z_g_XG.f(n_10, pjoint_Z_g_10, reZmat)$genZmat
    Z_01 = gen_Z_g_XG.f(n_01, pjoint_Z_g_01, reZmat)$genZmat
    Z_00 = gen_Z_g_XG.f(n_00, pjoint_Z_g_00, reZmat)$genZmat
    
    reXGmat = rbind(XGmat[xgInd_11==1,], XGmat[xgInd_10==1,],
                    XGmat[xgInd_01==1,], XGmat[xgInd_00==1,])
    trueZmat = rbind(Z_11, Z_10, Z_01, Z_00) 
    
    # generate T
    #system.time(trueTmat <- gen_T.f(K1, K2, K3, nsample, alp_1_vec, alp_2_vec, alp_3_vec, brks_1_vec, brks_2_vec, brks_3_vec, cor_mat))
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
    
    indata_s1 <- do.call("rbind", indata_s1.list)
    
    df_s2.list <- create_df_s2_jkjpkp.f(indata_s1.list, K1, K2, K3)
    
    est.prob.g <- mean(indata_s1$g)
    est.prob.x <- mean(indata_s1$x)
    
    est.pi_marg <- unlist(lapply(indata_s1.list, function(dt){
      mean(dt$truez)
    }))
    
    est.or <- unlist(lapply(df_s2.list, function(dt){
      
      est.pi11 <- mean(dt$truez_j*dt$truez_jp)
      est.pi10 <- mean(dt$truez_j*(1-dt$truez_jp))
      est.pi01 <- mean((1-dt$truez_j)*dt$truez_jp)
      est.pi00 <- mean((1-dt$truez_j)*(1-dt$truez_jp))
      
      est.or11 <- est.pi11*est.pi00/(est.pi10*est.pi01)
      
      est.or11
    }))
    
    
    est.survP <- unlist(lapply(indata_s1.list, function(dt){
      
      est.surv <- sum((dt$truet <= 1)*dt$truez)/sum(dt$truez)
      
    }))
    
    t.or <- c(or_11, or_22, or_33, or_12, or_13, or_23)
    t.survP <- c(rep(SurvP.ctr1, K1), rep(SurvP.ctr2, K2), rep(SurvP.ctr3, K3))
    t.pi_marg <- c(rep(pi1, K1), rep(pi2, K2), rep(pi3, K3))
    
    res <- (c(est.prob.g, est.prob.x, est.pi_marg, est.or, est.survP)- 
              c(prob.g, prob.x, t.pi_marg, t.or, t.survP))
    
    res.mat[i,] <- res
    
    print(i)
    
  }
  
  return(res.mat)
}

# test: succeed! (2021/01/18)
#res <- dataGen_check(nsim=100, nsample=2000, prob.g, prob.x, K1, K2, K3,  pi1, pi2, pi3, eta_11_vec, eta_21_vec, eta_31_vec,
#              SurvP.ctr1, SurvP.ctr2, SurvP.ctr3, #alp_11_vec, alp_21_vec, alp_31_vec,
#              kappa_1=1, kappa_2=1, kappa_3=1, or_11, or_22, or_33, or_12, or_13, or_23, 
#              tau_11, tau_22, tau_33, tau_12, tau_13, tau_23, visit.rate)
#colMeans(res)





