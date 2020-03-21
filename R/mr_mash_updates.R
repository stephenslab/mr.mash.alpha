###Update variational parameters, expected residuals, and ELBO components with or without scaling X
inner_loop_general <- function(X, rbar, mu, V, Vinv, w0, S0, ###note: V is only needed when not scaling X
                               precomp_quants, standardize, compute_ELBO, update_V){
  ###Create variables to store quantities
  n <- nrow(rbar)
  R <- ncol(rbar)
  p <- ncol(X)
  K <- length(S0)
  mu1   <- mu
  S1    <- array(0, c(R, R, p))
  w1    <- matrix(0, nrow=p, ncol=K)
  
  if(compute_ELBO){
    ##Initialize ELBO parameters
    var_part_tr_wERSS <- 0
    neg_KL <- 0
  }
  
  if(update_V){
    ##Initialize V parameters
    var_part_ERSS <- matrix(0, nrow=R, ncol=R)
  }
  
  ##Loop through the variables
  for(j in 1:p){
    
    #Remove j-th effect from expected residuals 
    rbar_j <- rbar + outer(X[, j], mu1[j, ])
    
    #Run Bayesian SLR
    if(standardize){
      bfit <- bayes_mvr_mix_scaled_X(X[, j], rbar_j, w0, S0, precomp_quants$S, precomp_quants$S1, 
                                     precomp_quants$SplusS0_chol, precomp_quants$S_chol)      
    } else {
      bfit <- bayes_mvr_mix_centered_X(X[, j], rbar_j, V, w0, S0, precomp_quants$xtx[j], Vinv, 
                                          precomp_quants$V_chol, precomp_quants$d, 
                                          precomp_quants$QtimesV_chol)
    }
    
    #Update variational parameters
    mu1[j, ]         <- bfit$mu1
    S1[, , j]        <- bfit$S1
    w1[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(compute_ELBO){
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
  if(compute_ELBO && update_V){
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL, var_part_ERSS=var_part_ERSS))
  } else if(compute_ELBO && !update_V){
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else if(!compute_ELBO && update_V) {
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1, var_part_ERSS=var_part_ERSS))
  } else { 
    return(list(rbar=rbar, mu1=mu1, S1=S1, w1=w1))
  }
}


###Perform one iteration of the outer loop with or without scaling X
mr_mash_update_general <- function(X, Y, mu1_t, V, Vinv, ldetV, w0, S0,
                                   precomp_quants, compute_ELBO, standardize, 
                                   update_V, version){
  ##Compute expected residuals
  rbar <- Y - X%*%mu1_t
  
  ##Update variational parameters, expected residuals, and ELBO components
  if(version=="R"){
    updates <- inner_loop_general(X=X, rbar=rbar, mu=mu1_t, V=V, Vinv=Vinv, w0=w0, S0=S0, 
                                  precomp_quants=precomp_quants, standardize=standardize,
                                  compute_ELBO=compute_ELBO, update_V=update_V)   
  } else if(version=="Rcpp"){
    updates <- inner_loop_general_rcpp_wrapper(X=X, Rbar=rbar, mu1=mu1_t, V=V, Vinv=Vinv, w0=w0,
                                              S0=simplify2array_custom(S0), precomp_quants=precomp_quants,
                                              standardize=standardize, compute_ELBO=compute_ELBO,
                                              update_V=update_V)
  }

  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  rbar    <- updates$rbar
  
  if(compute_ELBO && update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    var_part_ERSS <- updates$var_part_ERSS
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO, var_part_ERSS=var_part_ERSS))
  } else if(compute_ELBO && !update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, ELBO=ELBO))
  } else if(!compute_ELBO && update_V){
    var_part_ERSS <- updates$var_part_ERSS
    
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t, var_part_ERSS=var_part_ERSS))
  } else {
    return(list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t))
  }
}

###Wrapper for the Rcpp function to update variational parameters, expected residuals, 
###and ELBO components with or without scaling X
#' @importFrom Rcpp evalCpp
#' @useDynLib mr.mash.alpha
#' 
inner_loop_general_rcpp_wrapper <- function(X, Rbar, mu1, V, Vinv, w0, S0, precomp_quants, 
                                            standardize, compute_ELBO, update_V){

    out <- inner_loop_general_rcpp(X, Rbar, mu1, V, Vinv, w0, S0, precomp_quants, 
                                 standardize, compute_ELBO, update_V)

    ###Return output
  if(compute_ELBO && update_V){
    return(list(rbar=out$rbar, mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_tr_wERSS=out$var_part_tr_wERSS, 
                neg_KL=out$neg_KL, var_part_ERSS=out$var_part_ERSS))
  } else if(compute_ELBO && !update_V){
    return(list(rbar=out$rbar, mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_tr_wERSS=out$var_part_tr_wERSS, 
                neg_KL=out$neg_KL))
  } else if(!compute_ELBO && update_V) {
    return(list(rbar=out$rbar, mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_ERSS=out$var_part_ERSS))
  } else { 
    return(list(rbar=out$rbar, mu1=out$mu1, S1=out$S1, w1=out$w1))
  }
}