###Update variational parameters, expected residuals, and ELBO components
inner_loop <- function(X, rbar, mu, V, Vinv, w0, S0, xtx, V_chol, U0, d, Q){
  ###Create variables to store quantities
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1   <- mu
  S1    <- array(0, c(R, R, p))
  w1    <- matrix(0, nrow=p, ncol=K)
  
  if(!is.null(Vinv)){
    ##Initialize ELBO parameters
    var_part_tr_wERSS <- 0
    neg_KL <- 0
  }
  
  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1[j, ])
    
    #Run Bayesian SLR
    bfit <- bayes_mvr_mix_centered_X(X[, j], rbar_j, V, w0, S0, xtx[j], V_chol, U0, d, Q)
    
    #Update variational parameters
    mu1[j, ]         <- bfit$mu1
    S1[, , j]        <- bfit$S1
    w1[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(!is.null(Vinv)){
      var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%bfit$S1)*xtx[j])
      mu1_mat <- matrix(bfit$mu1, ncol=1)
      neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(X[, j], ncol=1), mu1_mat))+
                                             tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*(xtx[j])))
    }
    
    #Update expected residuals
    rbar <- rbar_j - outer(X[, j], mu1[j, ])
  }
  
  ###Return output
  if(!is.null(Vinv)){
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else {
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1))
  }
}

###Perform one iteration of the outer loop
mr_mash_update <- function(Y, X, mu1_t, w1_t, V, Vinv, ldetV, w0, S0, 
                           xtx, V_chol, U0, d, Q, update_w0, compute_ELBO){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  #Update w0 if requested
  if(update_w0 && !is.null(w1_t)){
    w0 <- update_weights(w1_t)
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  updates <- inner_loop(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=Vinv, w0=w0, S0=S0, xtx=xtx, V_chol=V_chol, U0=U0, d=d, Q=Q) 
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
}

###Update variational parameters, expected residuals, and ELBO components with scaled X
inner_loop_scaled_X <- function(X, rbar, mu, Vinv, w0, S0, S, S1, SplusS0_chol, S_chol, ldetSplusS0_chol, ldetS_chol){
  ###Create variables to store quantities
  n <- nrow(rbar)
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1c   <- mu
  S1c    <- array(0, c(R, R, p))
  w1c    <- matrix(0, nrow=p, ncol=K)
  
  if(!is.null(Vinv)){
    ##Initialize ELBO parameters
    var_part_tr_wERSS <- 0
    neg_KL <- 0
  }
  
  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1c[j, ])
    
    #Run Bayesian SLR
    bfit <- bayes_mvr_mix_scaled_X(X[, j], rbar_j, w0, S0, S, S1, SplusS0_chol, S_chol, ldetSplusS0_chol, ldetS_chol)
    
    #Update variational parameters
    mu1c[j, ]         <- bfit$mu1
    S1c[, , j]        <- bfit$S1
    w1c[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(!is.null(Vinv)){
      xtx <- n-1
      var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%bfit$S1)*xtx)
      mu1_mat <- matrix(bfit$mu1, ncol=1)
      neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(X[, j], ncol=1), mu1_mat))+
                                             tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*(xtx)))
    }
    
    #Update expected residuals
    rbar <- rbar_j - outer(X[, j], mu1c[j, ])
  }
  
  ###Return output
  if(!is.null(Vinv)){
    return(list(rbar=rbar, mu1=mu1c, S1=S1c, w1=w1c, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else {
    return(list(rbar=rbar, mu1=mu1c, S1=S1c, w1=w1c))
  }
}

###Perform one iteration of the outer loop with scaled X
mr_mash_update_scaled_X <- function(Y, X, mu1_t, w1_t, V, Vinv, ldetV, w0, S0, S, S1, 
                                    SplusS0_chol, S_chol, ldetSplusS0_chol, ldetS_chol, 
                                    update_w0, compute_ELBO){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  #Update w0 if requested
  if(update_w0 && !is.null(w1_t)){
    w0 <- update_weights(w1_t)
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  updates <- inner_loop_scaled_X(X=X, rbar=rbar, mu=mu1_t, Vinv=Vinv, w0=w0, S0=S0, 
                                 S=S, S1=S1, SplusS0_chol=SplusS0_chol, S_chol=S_chol, 
                                 ldetSplusS0_chol=ldetSplusS0_chol, ldetS_chol=ldetS_chol) 

  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
}


###Update variational parameters, expected residuals, and ELBO components with or without scaling X
inner_loop_general <- function(X, rbar, mu, V, Vinv, w0, S0, ###note: V is only needed when not scaling X
                               precomp_quants, standardize, update_V){
  ###Create variables to store quantities
  n <- nrow(rbar)
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1   <- mu
  S1    <- array(0, c(R, R, p))
  w1    <- matrix(0, nrow=p, ncol=K)
  
  if(!is.null(Vinv)){
    ##Initialize ELBO parameters
    var_part_tr_wERSS <- 0
    neg_KL <- 0
  }
  
  if(update_V){
    ##Initialize V parameters
    var_part_ERSS <- 0
  }
  
  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1[j, ])
    
    #Run Bayesian SLR
    if(standardize){
      bfit <- bayes_mvr_mix_scaled_X(X[, j], rbar_j, w0, S0, precomp_quants$S, precomp_quants$S1, 
                                     precomp_quants$SplusS0_chol, precomp_quants$S_chol, 
                                     precomp_quants$ldetSplusS0_chol, precomp_quants$ldetS_chol)      
    } else {
      bfit <- bayes_mvr_mix_centered_X(X[, j], rbar_j, V, w0, S0, precomp_quants$xtx[j], 
                                          precomp_quants$V_chol, precomp_quants$U0, precomp_quants$d, 
                                          precomp_quants$Q)
    }
    
    #Update variational parameters
    mu1[j, ]         <- bfit$mu1
    S1[, , j]        <- bfit$S1
    w1[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(!is.null(Vinv)){
      if(standardize){
        xtx <- n-1
      } else {
        xtx <- precomp_quants$xtx[j]
      }
      ELBO_parts <- compute_ELBO_terms(var_part_tr_wERSS, neg_KL, X[, j], rbar_j, bfit, xtx, Vinv)
      var_part_tr_wERSS <- ELBO_parts$var_part_tr_wERSS
      neg_KL <- ELBO_parts$neg_KL
    }
    
    #Compute V params
    if(update_V){
      if(standardize){
        xtx <- n-1
      } else {
        xtx <- precomp_quants$xtx[j]
      }
      var_part_ERSS <- compute_var_part_ERSS(var_part_ERSS, bfit, xtx)
    }
    
    #Update expected residuals
    rbar <- rbar_j - outer(X[, j], mu1[j, ])
  }
  
  ###Return output
  if(!is.null(Vinv) && update_V){
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL, var_part_ERSS=var_part_ERSS))
  } else if(!is.null(Vinv) && !update_V){
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else if(is.null(Vinv) && update_V) {
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_ERSS=var_part_ERSS))
  } else { 
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1))
  }
}


###Perform one iteration of the outer loop with or without scaling X
mr_mash_update_general <- function(X, Y, mu1_t, w1_t, V, Vinv, ldetV, w0, S0, var_part_ERSS,
                                   precomp_quants, update_w0, update_w0_method, 
                                   compute_ELBO, standardize, update_V){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  ##Update w0 if requested
  if(update_w0 && !is.null(w1_t)){
    if(update_w0_method=="EM"){
      w0 <- update_weights_em(w1_t)
    } else if(update_w0_method=="mixsqp"){
      w0em <- update_weights_em(w1_t)
      w0 <- update_weights_mixsqp(X=X, Y=Y, mu1_t=mu1_t, V=V, Vinv=Vinv, ldetV=ldetV, w0em=w0em, 
                                  S0=S0, precomp_quants=precomp_quants, standardize=standardize,
                                  stepsize.reduce = 0.5, stepsize.min = 1e-8)$w0
    }
  }
  
  ##Update V if requested
  if(update_V && !is.null(w1_t)){ #Use w1_t to check whether we are in the first iteration (if so, w1_t=NULL)
    V <- update_V(var_part_ERSS, rbar)
    if(standardize){
      precomp_quants <- precompute_quants_scaled_X(nrow(X), V, S0)
    } else {
      precomp_quants <- precompute_quants_centered_X(X, V, S0) 
    }
    if(compute_ELBO){
      Vinv <- chol2inv(precomp_quants$V_chol)
      ldetV <- chol2ldet(precomp_quants$V_chol)
    } 
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  updates <- inner_loop_general(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=Vinv, w0=w0, S0=S0, 
                                precomp_quants=precomp_quants, standardize=standardize,
                                update_V) 
  
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO && update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, V=V, Vinv=Vinv, ELBO=ELBO, var_part_ERSS=updates$var_part_ERSS))
  } else if(compute_ELBO && !update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else if(!compute_ELBO && update_V){
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, V=V, var_part_ERSS=updates$var_part_ERSS))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
}
