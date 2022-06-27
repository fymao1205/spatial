
# functions for the two-stage estimation 

# -----------------------------------------------------
# stage I: get MLE of alpj_mat, gamj_vec, gam_1, gam_2 
# -----------------------------------------------------

get.data.s1.f <- function(nsample, Jtype, KJjoint, XG, FT, Z, visit.rate){
  
  K_sum <- Jtype*KJjoint
  
  id <- rep(1:nsample, each=K_sum)
  x.vec <- rep(XG[,1], each=K_sum)
  g.vec <- rep(XG[,2], each=K_sum)
  
  type <- rep(rep(1:Jtype, each=KJjoint), nsample)
  joint <- rep(1:KJjoint, Jtype*nsample)
  truez <- Z
  ft <- FT
  
  # generate the last interval
  njoints <- length(ft)
  last.2visits <- gen.visit.f(visit.rate, njoints, ft)
  interval.u <- last.2visits[,1]
  interval.l <- last.2visits[,2]
  
  indata <- as.data.frame(cbind(id, type, joint, x.vec, g.vec, interval.l, interval.u,truez))
  colnames(indata) <- c("id", "type", "joint", "x", "g", "l", "r", "truez")
  
  return(indata)
  
}

Lcomp_1j <- function(l, r, alp, brks){
  res <- ppch(r, cuts=brks, levels=alp) - ppch(l, cuts=brks, levels=alp)
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

est.param.s1.f <- function(indata_s1, brks_list, gene.effect.type="common"){
  
  nsubj <- max(indata_s1$id)
  
  indata_type_1 <- indata_s1[indata_s1$type==1,]
  indata_type_2 <- indata_s1[indata_s1$type==2,]
  
  brks_1 <- brks_list[[1]]
  brks_2 <- brks_list[[2]]
  
  object.s1.f <- function(param){
    
    alp1_vec <- exp(param[1:2])
    alp2_vec <- exp(param[3:4])
    gam10 <- param[5]
    gam20 <- param[6]
    #gamj0_vec <- param[Jtype*Rj+1:Jtype]
    gam1 <- param[7]
    gam2 <- param[-(1:7)]
    
    if(gene.effect.type=="common"){
      
      pi_marg_1 <- expit.f(gam10+gam1*indata_type_1$x+gam2*indata_type_1$g)
      pi_marg_2 <- expit.f(gam20+gam1*indata_type_2$x+gam2*indata_type_2$g)
      
    }else if(gene.effect.type=="block"){
      
      pi_marg_1 <- expit.f(gam10+gam1*indata_type_1$x+gam2[1]*indata_type_1$g)
      pi_marg_2 <- expit.f(gam20+gam1*indata_type_2$x+gam2[2]*indata_type_2$g)
      
    }
    
    
    
    obs.logL_1 <- eval.obs.logL_j.f(indata_type_1$l, indata_type_1$r, alp1_vec, pi_marg_1, brks_1)
    obs.logL_2 <- eval.obs.logL_j.f(indata_type_2$l, indata_type_2$r, alp2_vec, pi_marg_2, brks_2)
    
    res <- -(sum(obs.logL_1)+sum(obs.logL_2))
    
    return(res)
    
  }
  
  if(gene.effect.type=="common"){
    dim.g <- 1
  }else if(gene.effect.type=="block"){
    dim.g <- 2
  }
  dim_p = 7+dim.g
  res.s1 <- nlm(f=object.s1.f, p=rep(0.05,dim_p), hessian=TRUE,gradtol=1e-06, steptol=1e-06)
  
  est.param <- res.s1$estimate
  est.alp1_vec <- exp(est.param[1:2])
  est.alp2_vec <- exp(est.param[3:4])
  est.gam10 <- est.param[5]
  est.gam20 <- est.param[6]
  est.gam1 <- est.param[7]
  est.gam2 <- est.param[-c(1:7)]
  s1 <- eval.obs.score_s1.f(indata_s1, c(est.alp1_vec, est.alp2_vec, 
                                         est.gam10, est.gam20, est.gam1, est.gam2),
                            brks_list, gene.effect.type=gene.effect.type)
  inv.A <- ginv(res.s1$hessian/nsubj)
  B <- t(s1) %*% s1/nsubj  #info_score$info/nsubj 
  inv.est.info <- (inv.A%*%B%*%inv.A)
  
  est.score <- colMeans(s1)  #res.s1$gradient
  est.ase <- sqrt(diag(inv.est.info)/nsubj)
  
  iter <- res.s1$iterations
  
  outdata_type_1 <- indata_type_1
  outdata_type_1$est.surv.l <- 1-ppch(outdata_type_1$l, cuts=brks_1, levels=est.alp1_vec)
  outdata_type_1$est.surv.r <- 1-ppch(outdata_type_1$r, cuts=brks_1, levels=est.alp1_vec)
  outdata_type_1$pi.marg <- expit.f(est.gam10+est.gam1*outdata_type_1$x+est.gam2*outdata_type_1$g)
  
  outdata_type_2 <- indata_type_2
  outdata_type_2$est.surv.l <- 1-ppch(outdata_type_2$l, cuts=brks_2, levels=est.alp2_vec)
  outdata_type_2$est.surv.r <- 1-ppch(outdata_type_2$r, cuts=brks_2, levels=est.alp2_vec)
  outdata_type_2$pi.marg <- expit.f(est.gam20+est.gam1*outdata_type_2$x+est.gam2*outdata_type_2$g)
  
  outdata <- as.data.frame(rbind(outdata_type_1, outdata_type_2))
  outdata <- outdata[order(outdata$id, outdata$type, outdata$joint),]
  outdata$inv.Phi.est.surv.l <- qnorm(outdata$est.surv.l)
  outdata$inv.Phi.est.surv.r <- qnorm(outdata$est.surv.r)
  
  return(list(res=res.s1,
              est.alp1_vec=est.alp1_vec,
              est.alp2_vec=est.alp2_vec,
              est.gam10=est.gam10,
              est.gam20=est.gam20,
              est.gam1=est.gam1,
              est.gam2=est.gam2,
              est.ase=est.ase,
              est.score=est.score,
              outdata=outdata,
              iter=iter
  ))	
  
}

# -----------------------------------------------------
# stage II: get MLE of sigma_jjp, pi_jjp
# -----------------------------------------------------

get.data.s2.f <- function(nsample, Jtype, KJjoint, outdata_s1){
  
  K_sum <- Jtype*KJjoint
  cm <- outdata_s1
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

# bivariate normal: CDF, PDF
bi.mvn.me <- function(inmat, sigma){
  
  mean <- rep(0,2)
  corr <- diag(2)
  corr[lower.tri(corr)] <- sigma
  corr[upper.tri(corr)] <- sigma
  
  if(dim(inmat)[1]==1){
    
    invec <- as.vector(inmat)
    if(sum(invec > .Machine$double.xmax)==1){
      return(pnorm(invec[!(invec > .Machine$double.xmax)]))
    }else{
      return(pmnorm(invec, varcov = corr))
    }
    
  }else{
    
    pos.Inf <- (inmat > .Machine$double.xmax)
    pos.oneInf.one.f <- (rowSums(pos.Inf)==1)
    if(any(pos.oneInf.one.f)){
      res <- rep(NA, dim(inmat)[1])
      inmat.y <- inmat[!pos.oneInf.one.f,]
      res.y <- pmnorm(inmat.y, varcov = corr)
      res[!pos.oneInf.one.f] <- res.y
      inmat.n <- inmat[pos.oneInf.one.f,]
      inmat.n[(inmat.n > .Machine$double.xmax)] <- 0
      if(sum(pos.oneInf.one.f)==1){
        invec.n <- sum(inmat.n)
      }else{
        invec.n <- as.vector(rowSums(inmat.n))
      }
      res.n <- pnorm(invec.n)
      res[pos.oneInf.one.f] <- res.n
    }else{
      res <- pmnorm(inmat, varcov = corr)
    }
    
    return( res)
  }
}

bi.dmvn.me.f <- function(inmat, sigma){
  
  mean <- rep(0,2)
  corr <- diag(2)
  corr[lower.tri(corr)] <- sigma
  corr[upper.tri(corr)] <- sigma
  
  res <- dmnorm(inmat, varcov = corr)
  res <- ifelse(is.nan(res),0,res)
  #res <- ifelse(is.infinite(res),0,res)
  return(as.double(res))
}

pi11_pw.f <- function(pi_marg_1, pi_marg_2, eta_12, cov){
  
  or_12 <- exp(cov %*% eta_12)
  if(!any(eta_12!=0)){
    return(pi_marg_1*pi_marg_2)
  }else{
    res <- 0.5*(or_12-1)^(-1)*(1+(pi_marg_1+pi_marg_2)*(or_12-1)-
                                 ((1+(pi_marg_1+pi_marg_2)*(or_12-1))^2
                                  +4*or_12*(1-or_12)*pi_marg_1*pi_marg_2)^(0.5))
  }
  return(res)
  #}
  
}

eval.obs.logL_jjp.f <- function(indata_s2_jjp, sigma_jjp, eta_jjp){
  
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
  
  #or_jjp <- exp(matrix(rep(1, dim(indata_s2_jjp)[1])) %*% eta_jjp)
  pi11_jjp <- pi11_pw.f(pi_marg_j, pi_marg_jp, eta_jjp, matrix(rep(1, dim(indata_s2_jjp)[1])))
  
  d2 <- bi.mvn.me(cbind(inv.l_j, inv.l_jp), sigma_jjp) - 
    bi.mvn.me(cbind(inv.l_j, inv.r_jp), sigma_jjp) - 
    bi.mvn.me(cbind(inv.r_j, inv.l_jp), sigma_jjp) +
    bi.mvn.me(cbind(inv.r_j, inv.r_jp), sigma_jjp)
  
  
  
  d1_j <- S.l_j - S.r_j
  d1_jp <- S.l_jp - S.r_jp
  
  pre_log <- d2*pi11_jjp + (1-D_jp)*d1_j*(pi_marg_j-pi11_jjp)+
    (1-D_j)*d1_jp*(pi_marg_jp-pi11_jjp)+ 
    (1-D_jp-D_j+D_j*D_jp)*(1-pi_marg_j-pi_marg_jp+pi11_jjp)
  
  #res <- (log(pre_log))
  res <- aggregate((log(pre_log)), list(indata_s2_jjp$id), sum, drop = T, simplify = T)[, -1]
  
  res <- ifelse(is.nan(res), -Inf, res)
  #res <- aggregate(res, list(indata_s2_jjp$id), sum, drop = T, simplify = T)[, -1]
  return(res)
  
}

est.param.s2.f <- function(Kj_vec, indata_s2, pair.asso.type=NULL, s2.num.cores=2){
  
  J <- length(Kj_vec)
  K_2 <- choose(sum(Kj_vec),2)
  if(is.null(pair.asso.type)){
    pair.asso.type <- "common"
  }
  
  #jjp0.mat <- unique(indata_s2[,c("type_j","joint_j","type_jp","joint_jp")])
  
  if(pair.asso.type=="common"){
    num_param <- 1
    object.f <- function(param){
      
      sig_vec <- atan(param[(1:num_param)])*2/pi
      eta_vec <- param[-((1:num_param))]
      
      t1 <- Sys.time()
      obs.logL_s2 <- eval.obs.logL_jjp.f(indata_s2, sig_vec, eta_vec)
      t2 <- Sys.time()
      print(t2-t1)
      
      res <- -sum(unlist(obs.logL_s2))
      #res <- ifelse(is.nan(res), Inf, res)
      return(res)
      
    }
    
    res.s2 <- nlm(f=object.f, p=rep(0.5,2*num_param),hessian=TRUE,gradtol=1e-05, steptol=1e-06)
    est.param <- res.s2$estimate  #as.vector(res.s2$par) #
    iter <- res.s2$iterations #res.s2$counts[1]#
    
    est.sigma_pw <- atan(est.param[(1:num_param)])*2/pi
    est.eta_pw <- est.param[-(1:num_param)]
    
  }else if(pair.asso.type=="jjp"){
    
    indata_s2$group <- as.factor(paste0(indata_s2$type_j, 
                                        #indata_s2$joint_j, 
                                        indata_s2$type_jp
                                        #,indata_s2$joint_jp
    ))
    indata_s2.list <- split(indata_s2, indata_s2$group)
    
    object.f <- function(param, indata_jjp){
      
      print(param)
      sig <- atan(param[1])*2/pi
      eta <- param[2]
      
      t1 <- Sys.time()
      indata_jjp$group <- as.factor(paste0(indata_jjp$joint_j,indata_jjp$joint_jp))
      indt_jjp.list <- split(indata_jjp, indata_jjp$group)
      res <-  mclapply(indt_jjp.list, eval.obs.logL_jjp.f,
                       sig=sig, eta=eta, mc.cores=2)
      
      t2 <- Sys.time()
      print(t2-t1)
      res <- -sum(unlist(res))
      print(res)
      return(res)
      
    }
    
    system.time(
      res.list <-lapply(indata_s2.list, function(dt){
        nlm(f=object.f, p=rep(0.05,2),indata_jjp=dt,hessian=TRUE,gradtol=1e-05, steptol=1e-06)
      }
      #, mc.cores=s2.num.cores
      )
    )
    
    est <- unlist(lapply(res.list, function(x){
      x$estimate
    }))
    hess <- matrix(0, nrow=2*length(res.list), ncol=2*length(res.list))
    for(j in 1:length(res.list)){
      hess[2*(j-1)+1:2,2*(j-1)+1:2] <- res.list[[j]]$hessian
    }
    
    grad <- unlist(lapply(res.list, function(x){
      x$grad
    }))
    
    iter <- max(unlist(lapply(res.list, function(x){
      x$iter
    })))
    
    res.s2 <- list(estimate=est, 
                   hessian=hess,
                   gradient=grad,
                   iterations=iter)
    
    est.param <- est  #as.vector(res.s2$par) #
    num_param <- length(res.list)
    est.sigma_pw <- atan(est.param[(1:num_param)*2-1])*2/pi
    est.eta_pw <- est.param[(1:num_param)*2]
    
  }else if(pair.asso.type=="full"){
    num_param <- K_2
    jjp0.mat <- unique(indata_s2[,c("type_j","joint_j","type_jp","joint_jp")])
    
    object.f <- function(param){
      
      sig_vec <- atan(param[(1:num_param)])*2/pi
      eta_vec <- param[-((1:num_param))]
      
      jjp.mat <- cbind(jjp0.mat, sig_vec, eta_vec)
      
      t1 <- Sys.time()
      obs.logL_s2 <- mclapply(1:dim(jjp.mat)[1], function(i){
        jjp <- as.double(jjp.mat[i,])
        indata_s2_jjp <- indata_s2.list[[i]]
        eval.obs.logL_jjp.f(indata_s2_jjp, jjp[dim(jjp0.mat)[2]+1], jjp[dim(jjp0.mat)[2]+2])
      }, mc.cores=s2.num.cores)
      
      t2 <- Sys.time()
      print(t2-t1)
      
      res <- -sum(unlist(obs.logL_s2))
      #res <- ifelse(is.nan(res), Inf, res)
      return(res)
      
    }
    
    res.s2 <- nlm(f=object.f, p=rep(0.5,2*num_param),hessian=TRUE,gradtol=1e-05, steptol=1e-06)
    est.param <- res.s2$estimate  #as.vector(res.s2$par) #
    iter <- res.s2$iterations #res.s2$counts[1]#
    
    est.sigma_pw <- atan(est.param[(1:num_param)])*2/pi
    est.eta_pw <- est.param[-(1:num_param)]
  }
  
  
  
  return(list(res=res.s2,
              est.sigma_pw=est.sigma_pw,
              est.eta_pw=est.eta_pw,
              iter=iter))
}

