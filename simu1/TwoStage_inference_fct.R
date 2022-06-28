
# score and information estimation for inference

eval.obs.score_s1.f <- function(indata_s1, theta1, brks_list, gene.effect.type="common"){
  
  if(is.null(gene.effect.type)){ gene.effect.type = "common"}
  
  Jtype <- length(brks_list)
  the1 <- theta1
  alp.list <-vector(mode="list", length=Jtype)
  
  for(j in 1:Jtype){
    len_j <- length(brks_list[[j]])+1
    alp.list[[j]] <- the1[1:len_j]
    the1 <- the1[-(1:len_j)]
  }
  
  gam_j0_vec <- the1[1:Jtype]
  if(gene.effect.type=="common"){
    gam_sup <- matrix(rep(the1[-(1:Jtype)], Jtype), nrow=Jtype, byrow=1)
  }else if(gene.effect.type=="block"){
    gam_gj_vec <- the1[(length(the1)-Jtype+1):length(the1)]
    gam_sup <- cbind(rep(the1[- c(1:Jtype, (length(the1)-Jtype+1):length(the1))], Jtype),
                     gam_gj_vec)
  }
  
  gam <- cbind(gam_j0_vec, gam_sup)
  
  res.list <- mclapply(1:Jtype, function(j){
    indata_type_j <- indata_s1[indata_s1$type==j,]
    brks_j <- brks_list[[j]]
    outdata_res_j <- update_indata_s1.f(indata_type_j, brks_j, alp.list[[j]], 
                                        gam[j,1], gam[j,2], gam[j,3])
    outdata_j <- outdata_res_j$outdata_j
    ext_outdata_j <- outdata_res_j$ext_outdata_j
    eval.obs.score_j <- eval.obs.score.s1_j.f(outdata_j, ext_outdata_j)
  })
  
  res <- do.call("cbind", res.list)
  
  Rj_vec <- unlist(lapply(1:length(brks_list), function(j){
    length(brks_list[[j]])+1
  }))
  
  start <-cumsum(Rj_vec)
  pos.gam.vec <- unlist(lapply(1:Jtype,function(j){
    start[[j]]+(j-1)*3 + 1:3
  }))
  
  if(gene.effect.type=="common"){
    res <- cbind(res[,-pos.gam.vec], res[,pos.gam.vec[c(1,4)]], 
                 res[, pos.gam.vec[c(2:3)]]+res[, pos.gam.vec[c(5:6)]])
  }else if(gene.effect.type=="block"){
    res <- cbind(res[,-pos.gam.vec], res[,pos.gam.vec[c(1,4)]], 
                 res[, pos.gam.vec[c(2)]]+res[, pos.gam.vec[c(5)]],
                 res[,pos.gam.vec[c(3,6)]])
  }
  
  #info <- t(eval.obs.score) %*% eval.obs.score
  #sqrt(diag(ginv(info/nsubj))/nsubj)
  
  #return(list(info=info,
  #            score=eval.obs.score))
  return(res)  
}

d_d2_d_sig.f <- function(sigma_12,inv.Phi.S.l_1, inv.Phi.S.r_1, 
                         inv.Phi.S.l_2, inv.Phi.S.r_2){
  
  tr.sig_12 <- tan(sigma_12*pi/2)
  res <- (bi.dmvn.me.f(cbind(inv.Phi.S.l_1, inv.Phi.S.l_2), sigma_12) - 
            bi.dmvn.me.f(cbind(inv.Phi.S.l_1, inv.Phi.S.r_2), sigma_12) - 
            bi.dmvn.me.f(cbind(inv.Phi.S.r_1, inv.Phi.S.l_2), sigma_12) +
            bi.dmvn.me.f(cbind(inv.Phi.S.r_1, inv.Phi.S.r_2), sigma_12))*2/pi*1/(1+tr.sig_12^2)
  
  return(as.double(res))
}

d_pi11_pw_d_eta_12.f <- function(pi_marg_1, pi_marg_2, eta_12, cov){
  
  or_12 <- as.vector(exp(cov %*% eta_12))
  if(eta_12==0){
    return(rep(0, length(eta_12)))
  }else{
  
    phi_12 <- ((1+(pi_marg_1+pi_marg_2)*(or_12-1))^2
               +4*or_12*(1-or_12)*pi_marg_1*pi_marg_2)^(0.5)
    pre.res <- or_12*(2*(or_12-1)^2*phi_12)^(-1)*(1-phi_12+(or_12-1)*(pi_marg_1+pi_marg_2-2*pi_marg_1*pi_marg_2))
    res.list <- lapply(1:dim(cov)[1], function(i){
      cov[i,] * pre.res[i]
    })
    #cov * pre.res
    res <- do.call("rbind", res.list)
    return(res)
  }
  return(as.vector(res))
}

eval.obs.score_s2.i12.f <- function(indata_s2_jjp, sigma_jjp, eta_jjp){
  
  D_j <- is.finite(indata_s2_jjp$r_j)
  D_jp <- is.finite(indata_s2_jjp$r_jp)
  S.l_j <- indata_s2_jjp$est.surv.l_j
  S.r_j <- indata_s2_jjp$est.surv.r_j
  pi_marg_j <- indata_s2_jjp$pi_marg_j
  inv.l_j <- indata_s2_jjp$inv.Phi.est.surv.l_j
  inv.r_j <- indata_s2_jjp$inv.Phi.est.surv.r_j
  S.l_jp <- indata_s2_jjp$est.surv.l_jp
  S.r_jp <- indata_s2_jjp$est.surv.r_jp
  pi_marg_jp <- indata_s2_jjp$pi_marg_jp
  inv.l_jp <- indata_s2_jjp$inv.Phi.est.surv.l_jp
  inv.r_jp <- indata_s2_jjp$inv.Phi.est.surv.r_jp
  
  pi11_jjp <- pi11_pw.f(pi_marg_j, pi_marg_jp, eta_jjp, matrix(rep(1, dim(indata_s2_jjp)[1])))
  
  d2 <- bi.mvn.me(cbind(inv.l_j, inv.l_jp), sigma_jjp) - 
    bi.mvn.me(cbind(inv.l_j, inv.r_jp), sigma_jjp) - 
    bi.mvn.me(cbind(inv.r_j, inv.l_jp), sigma_jjp) +
    bi.mvn.me(cbind(inv.r_j, inv.r_jp), sigma_jjp)
  
  d_d2_d_sig <- d_d2_d_sig.f(sigma_jjp,inv.l_j, inv.r_j, 
                             inv.l_jp, inv.r_jp)
  
  d_pi11_pw_d_eta_12 <- d_pi11_pw_d_eta_12.f(pi_marg_j, pi_marg_jp,
                                             eta_jjp, matrix(rep(1, dim(indata_s2_jjp)[1])))
  
  d1_j <- S.l_j - S.r_j
  d1_jp <- S.l_jp - S.r_jp
  
  pre_log <- d2*pi11_jjp + (1-D_jp)*d1_j*(pi_marg_j-pi11_jjp)+
    (1-D_j)*d1_jp*(pi_marg_jp-pi11_jjp)+ 
    (1-D_jp-D_j+D_j*D_jp)*(1-pi_marg_j-pi_marg_jp+pi11_jjp)
  
  res.sig12 <- d_d2_d_sig*pi11_jjp/pre_log
  res.eta12 <- (d2*d_pi11_pw_d_eta_12 + (1-D_jp)*d1_j*(-d_pi11_pw_d_eta_12)
                +(1-D_j)*d1_jp*(-d_pi11_pw_d_eta_12) +
                  (1-D_jp-D_j+D_j*D_jp)*(d_pi11_pw_d_eta_12))/pre_log
  
  res.sig12 <- aggregate(res.sig12, list(indata_s2_jjp$id), sum, drop = T, simplify = T)[, -1]
  res.eta12 <- aggregate(res.eta12, list(indata_s2_jjp$id), sum, drop = T, simplify = T)[, -1]
  res <- cbind(res.sig12, res.eta12)
  res <- ifelse(is.nan(res), -Inf, res)
  return(res)
  
}

eval.obs.score_s2.f <- function(indata_s2, sig_vec, eta_vec, pair.asso.type, ase.num.cores){
  
  indata_s2$group <- as.factor(paste0(indata_s2$type_j, #indata_s2$joint_j, 
                                      indata_s2$type_jp #, indata_s2$joint_jp
  ))
  indata_s2.list <- split(indata_s2, indata_s2$group)
  
  if(pair.asso.type=="full"){
    jjp0.mat <- unique(indata_s2[,c("type_j","joint_j","type_jp","joint_jp")])
    jjp.mat <- cbind(jjp0.mat, sig_vec, eta_vec)
    
    res.list <- mclapply(1:dim(jjp.mat)[1], function(i){
      jjp <- as.double(jjp.mat[i,])
      indata_s2_jjp <- indata_s2.list[[i]]
      eval.obs.score_s2.i12.f(indata_s2_jjp, jjp[5], jjp[6])
    }, mc.cores=ase.num.cores)
    
    res <- do.call("cbind", res.list)
    
  }else if(pair.asso.type=="jjp"){
    
    jjp0.mat <- unique(indata_s2[,c("type_j","type_jp")])
    jjp.mat <- cbind(jjp0.mat, sig_vec, eta_vec)
    
    res.list <- mclapply(1:dim(jjp.mat)[1], function(i){
      jjp <- as.double(jjp.mat[i,])
      indata_s2_jjp <- indata_s2.list[[i]]
      eval.obs.score_s2.i12.f(indata_s2_jjp, jjp[3], jjp[4])
    }, mc.cores=ase.num.cores)
    res <- do.call("cbind", res.list)
  }else{
    res <- eval.obs.score_s2.i12.f(indata_s2, sig_vec, eta_vec)
  }
  
  #res <- res[,c((1:(dim(res)[2]/2))*2-1, (1:(dim(res)[2]/2))*2)] #c(res[(1:(length(res)/2))*2-1], res[-((1:(length(res)/2))*2-1)])
  return(res)
  
}


get_indata_s2.f <- function(Kj_vec, theta1, brks.list, indata){
  
  nsample <- length(unique(indata$id))
  
  J <- length(brks.list)
  piece.sum <- length(unlist(brks.list)) + J
  alpha.list <- list()
  from <- 0 
  for(j in 1:J){
    to <- from + length(brks.list[[j]])+1
    alpha.list[[j]] <- theta1[(from+1):to]
    from <- to 
  }
  
  gam <- matrix(NA, nrow=J, ncol=3)
  gam[,1] <- theta1[(piece.sum)+1:J]
  gam[,2] <- rep(theta1[(piece.sum+J)+1],J)
  #if(gene.effect.type=="common"){
  #  gam[,3] <- rep(theta1[(piece.sum+J)+2],J)
  #}else if(gene.effect.type=="block"){
    gam[,3] <- theta1[-(1:((piece.sum+J)+1))]
  #}
  
  
  indata_type_j.list <- lapply(1:length(Kj_vec), function(j){
    indata_type_j <- indata[indata$type==j,]
    brks_j <- brks.list[[j]]
    outdata_type_j <- indata_type_j
    outdata_type_j$est.surv.l <- 1-ppch(outdata_type_j$l, cuts=brks.list[[j]], levels=alpha.list[[j]])
    outdata_type_j$est.surv.r <- 1-ppch(outdata_type_j$r, cuts=brks.list[[j]], levels=alpha.list[[j]])
    outdata_type_j$pi.marg <- expit.f(gam[j,1]+gam[j,2]*outdata_type_j$x+gam[j,3]*outdata_type_j$g)
    outdata_type_j
  })
  
  outdata <- as.data.frame(do.call("rbind", indata_type_j.list))
  outdata <- outdata[order(outdata$id, outdata$type, outdata$joint),]
  outdata$inv.Phi.est.surv.l <- qnorm(outdata$est.surv.l)
  outdata$inv.Phi.est.surv.r <- qnorm(outdata$est.surv.r)
  
  K_sum <- sum(Kj_vec)
  cm <- outdata
  cm.ext <- cm[rep(1:nrow(cm), times=rep((K_sum-1:K_sum),nsample)),]
  cmp <- cm[,-c(4,5)]
  cmp.ext <-  do.call(rbind, lapply(1:nsample, function(ith){
    cmp_ith <- cmp[cmp$id==ith,]
    do.call(rbind,rep(list(cmp_ith), K_sum))
  }))
  indexp <- rep(as.vector(unlist(sapply(1:K_sum,function(x){
    -(1:x)-K_sum*(x-1)
  }))), times=nsample) - (rep(1:nsample,each=sum(1:K_sum))-1)*K_sum^2
  cmp.ext <- cmp.ext[indexp,]
  data_s2.mat <- cbind(cm.ext, cmp.ext[,-c(1)])
  
  indata_s2 <- as.data.frame(data_s2.mat)
  colnames(indata_s2) <- c("id", "type_j", "joint_j", "x", "g","l_j", "r_j", "truez_j", "est.surv.l_j", "est.surv.r_j", "pi_marg_j", "inv.Phi.est.surv.l_j", "inv.Phi.est.surv.r_j", 
                           "type_jp", "joint_jp", "l_jp", "r_jp", "truez_jp", "est.surv.l_jp", "est.surv.r_jp", "pi_marg_jp", "inv.Phi.est.surv.l_jp", "inv.Phi.est.surv.r_jp")
  
  return(indata_s2)
  
}

eval.A11.f <- function(delta, indata_s1, alpha.list, brks.list, gam, gene.effect.type="common",
                       ase.num.cores){
  
  theta1 <- c(unlist(alpha.list), gam[,1], gam[1,-1])
  len <- length(theta1)
  
  dU1.list <- mclapply(1:len, function(l){
    the1.p <- the1.n <- theta1
    the1.p[l] <- theta1[l]+delta
    the1.n[l] <- theta1[l]-delta
    alp.list.p <- 
    tmpU1.p <- eval.obs.score_s1.f(indata_s1, the1.p, brks_list,gene.effect.type=gene.effect.type)
    tmpU1.n <- eval.obs.score_s1.f(indata_s1, the1.n, brks_list,gene.effect.type=gene.effect.type)
    tmpdU1 <- (tmpU1.p - tmpU1.n)/(2*delta)
    colSums(tmpdU1)
  },mc.cores=ase.num.cores)
  
  dU1 <- do.call("cbind", dU1.list)
  return(dU1)
}

eval.A12.f <- function(delta, Kj_vec, indata, alpha.list, brks.list, gam, gene.effect.type="common",
                       sig_vec, eta_vec, pair.asso.type,ase.num.cores){
  
  K_2 <- combn(rep(1:length(Kj_vec),Kj_vec), 2)
  
  if(gene.effect.type=="common"){
    gam_g <- unique(gam[,3])
  }else if(gene.effect.type=="block"){
    gam_g <- gam[,3]
  }
  theta1 <- c(unlist(alpha.list), 
              gam[,1], gam[1,2], gam_g)
  len <- length(theta1)
  
  dU2.list <- mclapply(1:len, function(l){
    the1.p <- the1.n <- theta1
    the1.p[l] <- theta1[l]+delta
    the1.n[l] <- theta1[l]-delta
    indata_s2.p <- get_indata_s2.f(Kj_vec, the1.p, brks.list, indata)
    indata_s2.n <- get_indata_s2.f(Kj_vec, the1.n, brks.list, indata)
    tmpU2.p <- eval.obs.score_s2.f(indata_s2.p, sig_vec, eta_vec, pair.asso.type,ase.num.cores)
    tmpU2.n <- eval.obs.score_s2.f(indata_s2.n, sig_vec, eta_vec, pair.asso.type,ase.num.cores)
    tmpdU2 <- (tmpU2.p - tmpU2.n)/(2*delta)
    colSums(tmpdU2)
  },mc.cores=ase.num.cores)
  
  dU2 <- do.call("cbind", dU2.list)
  return(dU2)
}

eval.B_ts.f <- function(indata_s1, alp.list, brks_list, gam, gene.effect.type="common",
                        indata_s2, sig_vec, eta_vec, pair.asso.type,ase.num.cores){
  
  theta1 <- c(unlist(alp.list), gam[,1], gam[1,-1])
  s1 <- eval.obs.score_s1.f(indata_s1, theta1, brks_list,gene.effect.type)
  s2 <- eval.obs.score_s2.f(indata_s2, sig_vec, eta_vec, pair.asso.type,ase.num.cores)
  s <- cbind(s1,s2)
  #res.mat <- apply(s, 1, function(x){(x) %*% t(x)})
  #res <- matrix(rowMeans(res.mat), nrow=dim(s)[2])
  res <- t(s) %*% s
  return(res)
}

eval.asvar_ts.f <- function(delta, nsample, Kj_vec, est.A11, est.A22,
                            indata_s1, est.alpha.list, brks_list, 
                            est.gam, gene.effect.type="common",
                            indata_s2, est.sig_vec, est.eta_vec, pair.asso.type,ase.num.cores){
  
  # A11
  # A12
  system.time(
    A12 <- eval.A12.f(delta=delta, Kj_vec=Kj_vec, 
                      indata=indata_s1, alpha.list=est.alpha.list, brks.list=brks_list, 
                      gam=est.gam, gene.effect.type=gene.effect.type,
                      sig_vec=est.sig_vec, eta_vec=est.eta_vec, pair.asso.type=pair.asso.type,
                      ase.num.cores=ase.num.cores)
  )
  # A22
  if(gene.effect.type=="common"){
    dim.g <- 1
  }else if(gene.effect.type=="block"){
    dim.g <- length(brks_list)
  }
  len1 <- length(c(unlist(est.alpha.list), est.gam[,1]))+dim.g+1
  len <- len1 + length(c(est.sig_vec, est.eta_vec))
  A <- matrix(0, len, len)
  
  A[1:len1, 1:len1] <- est.A11
  A[(len1+1):len,1:len1] <- A12/nsample
  A[(len1+1):len,(len1+1):len] <- est.A22
  
  B <- eval.B_ts.f(indata_s1, est.alpha.list, brks_list, est.gam,gene.effect.type=gene.effect.type,
                   indata_s2, est.sig_vec, est.eta_vec, pair.asso.type,ase.num.cores)
  
  inv.A <- try(ginv(A/nsample), silent=TRUE)
  if(class(inv.A) == "try-error"){
    warning("Fail to inverse of information matrix")
    V <- NA
  }else{
    
    B <- B/nsample
    V <- inv.A %*% B %*% (inv.A)/nsample
    
  }
  
  return(V)
  
}

eval.obs.score.s1_ijk.f <- function(outdata_jk, ext_outdata_jk){
  
  delta_jk <- is.finite(outdata_jk$r)*1
  denom <- ((ext_outdata_jk$surv_l - ext_outdata_jk$surv_r)*ext_outdata_jk$pi_j +(1-ext_outdata_jk$pi_j))
  score_alp_j.mat <- matrix(0, nrow=length(delta_jk), ncol=max(ext_outdata_jk$piece_j))
  score_alp_j.mat[delta_jk==1,] <- matrix(((-ext_outdata_jk$surv_l*ext_outdata_jk$wl*ext_outdata_jk$alp_j
                                            +ext_outdata_jk$surv_r*ext_outdata_jk$wr*ext_outdata_jk$alp_j)/(ext_outdata_jk$surv_l - ext_outdata_jk$surv_r)), nrow=length(delta_jk), byrow=1)[delta_jk==1,]
  score_alp_j.mat[delta_jk==0,] <- matrix(((-ext_outdata_jk$surv_l*ext_outdata_jk$wl*ext_outdata_jk$alp_j+ 0 )*ext_outdata_jk$pi_j/denom), nrow=length(delta_jk), byrow=1)[delta_jk==0,]
  
  score_gam_j.mat <- matrix(0, nrow=length(delta_jk), ncol=3)
  cov.mat <- cbind(rep(1,length(delta_jk)), outdata_jk$x, outdata_jk$g)#matrix(rep(c(1,x,g), length(delta_jk)), ncol=length(alp_j), byrow=1)
  score_gam_j.mat[delta_jk==1,] <- ((1-outdata_jk$pi_j)*cov.mat)[delta_jk==1,]
  
  pre_score_j <- (outdata_jk$surv_l-outdata_jk$surv_r-1)*outdata_jk$pi_j*(1-outdata_jk$pi_j)/((outdata_jk$surv_l - outdata_jk$surv_r)*outdata_jk$pi_j +(1-outdata_jk$pi_j)) #as.vector(matrix((surv_ljk-surv_rjk-1)*pi_j*(1-pi_j)/denom, nrow=length(delta_jk), byrow=1)[,1])
  score_gam_j.mat[delta_jk==0,] <- (pre_score_j*cov.mat)[delta_jk==0,]
  
  score_mat <- cbind(score_alp_j.mat, score_gam_j.mat)
  
  return(score_mat)
}

eval.obs.score.s1_j.f <- function(outdata_j, ext_outdata_j){
  
  R_j <- dim(ext_outdata_j)[1]/dim(outdata_j)[1]
  K_j <- max(outdata_j$joint)
  
  res.list <- mclapply(1:K_j, function(k){
    ext_outdata_jk <- ext_outdata_j[ext_outdata_j$joint==k,]
    outdata_jk <- outdata_j[outdata_j$joint==k,]
    res <- eval.obs.score.s1_ijk.f(outdata_jk, ext_outdata_jk)
    return(res)
  }, mc.cores=2)
  
  res <- matrix(0, nrow=dim(res.list[[1]])[1], ncol=dim(res.list[[1]])[2])
  for(k in 1:K_j){
    res <- res + res.list[[k]]
  }
  
  #  res1 <- (outdata_j$surv_l - outdata_j$surv_r)*outdata_j$pi_j + (1-outdata_j$pi_j)*(1-is.finite(outdata_j$r))
  #  z2 <- -(ext_outdata_j$wr*ext_outdata_j$alp_j)*ext_outdata_j$surv_r 
  #  z2 <- ifelse(is.finite(ext_outdata_j$r), z2, 0)
  #  res2_vec.mat <- cbind(matrix(ext_outdata_j$pi_j*(-(ext_outdata_j$wl*ext_outdata_j$alp_j)*ext_outdata_j$surv_l - z2), ncol = R_j, byrow = 1), 
  #                ((outdata_j$surv_l - outdata_j$surv_r)-(1-is.finite(outdata_j$r)))*outdata_j$pi_j*(1-outdata_j$pi_j)*cbind(rep(1,dim(outdata_j)[1]),outdata_j$x,outdata_j$g) )
  #  res <- matrix(as.vector(res2_vec.mat), nrow = length(res1))/res1
  
  return(res)
}

update_indata_s1.f <- function(indata_type_j, brks_j, alp_vec_j, 
                               gam_j0, gam_1, gam_2){
  
  R_j <- (length(brks_j)+1)
  
  indata_type_j$pi_j <- expit.f(gam_j0 + gam_1*indata_type_j$x + gam_2*indata_type_j$g)
  indata_type_j$surv_l <- 1-ppch(indata_type_j$l, cuts=brks_j, levels=alp_vec_j)
  indata_type_j$surv_r <- 1-ppch(indata_type_j$r, cuts=brks_j, levels=alp_vec_j)
  
  ext_indata_j <- t(apply(indata_type_j, 1, function(x){
    replicate(R_j, x)
  }))
  
  ext_outdata_j <- as.data.frame(rbind(ext_indata_j[,1:dim(indata_type_j)[2]],
                                       ext_indata_j[,1:dim(indata_type_j)[2]+dim(indata_type_j)[2]]))
  colnames(ext_outdata_j) <- colnames(indata_type_j)
  ext_outdata_j <- ext_outdata_j[order(ext_outdata_j$id, ext_outdata_j$type, ext_outdata_j$joint),]
  ext_outdata_j$piece_j <- rep(1:R_j, dim(indata_type_j)[1])
  ext_outdata_j$alp_j <- rep(alp_vec_j, dim(indata_type_j)[1])
  ext_outdata_j$gam_j0 <- rep(gam_j0, dim(ext_outdata_j)[1])
  ext_outdata_j$wl <- as.vector(unlist(mclapply(indata_type_j$l, function(l){
    c(l, brks_j)[order(c(l, brks_j))] - c(0, brks_j)
  }, mc.cores=2)))
  ext_outdata_j$wr <- as.vector(unlist(mclapply(indata_type_j$r, function(u){
    c(u, brks_j)[order(c(u, brks_j))] - c(0, brks_j)
  }, mc.cores=2)))
  
  res <- list(outdata_j=indata_type_j, 
              ext_outdata_j = ext_outdata_j)
  return(res)
}





