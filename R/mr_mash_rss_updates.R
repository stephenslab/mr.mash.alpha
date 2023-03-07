###Update variational parameters, expected residuals, and ELBO components with or without scaling X
inner_loop_general_rss_R <- function(n, XtY, XtXmu1, mu1, V, Vinv, w0, S0, ###note: V is only needed when not scaling X
                               precomp_quants, standardize, compute_ELBO, update_V,
                               update_order, eps){
  ###Create variables to store quantities
  r <- ncol(mu1)
  p <- nrow(mu1)
  K <- length(S0)
  S1    <- array(0, c(r, r, p))
  w1    <- matrix(0, nrow=p, ncol=K)
  var_part_tr_wERSS <- 0
  neg_KL <- 0
  var_part_ERSS <- matrix(0, nrow=r, ncol=r)
  
  ##Loop through the variables
  for(j in update_order){
    
    if(standardize){
      xtx <- n-1
    } else {
      xtx <- precomp_quants$xtx[j]
    }
    
    #Remove j-th effect from expected residuals 
    xtRbar_j <- XtY[j, ] - XtXmu1[j, ] + xtx*mu1[j, ]
    
    #Run Bayesian SLR
    if(standardize){
      bfit <- bayes_mvr_mix_standardized_X_rss(n, xtRbar_j, w0, S0, precomp_quants$S, precomp_quants$S1, 
                                     precomp_quants$SplusS0_chol, precomp_quants$S_chol, eps)      
    } else {
      bfit <- bayes_mvr_mix_centered_X_rss(xtRbar_j, V, w0, S0, xtx, Vinv, 
                                          precomp_quants$V_chol, precomp_quants$d, 
                                          precomp_quants$QtimesV_chol, eps)
    }
    
    #Update variational parameters
    mu1[j, ]         <- bfit$mu1
    S1[, , j]        <- bfit$S1
    w1[j, ]          <- bfit$w1
    
    #Compute ELBO params
    if(compute_ELBO){
      ELBO_parts <- compute_ELBO_rss_terms(var_part_tr_wERSS, neg_KL, xtRbar_j, bfit, xtx, Vinv)
      var_part_tr_wERSS <- ELBO_parts$var_part_tr_wERSS
      neg_KL <- ELBO_parts$neg_KL
    }
    
    #Compute V params
    if(update_V){
      var_part_ERSS <- compute_var_part_ERSS(var_part_ERSS, bfit, xtx)
    }
  }
  
  ###Return output
  if(compute_ELBO && update_V){
    return(list(mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL, var_part_ERSS=var_part_ERSS))
  } else if(compute_ELBO && !update_V){
    return(list(mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else if(!compute_ELBO && update_V) {
    return(list(mu1=mu1, S1=S1, w1=w1, var_part_ERSS=var_part_ERSS))
  } else { 
    return(list(mu1=mu1, S1=S1, w1=w1))
  }
}

### Wrapper for the Rcpp function to update variational parameters,
### expected residuals, and ELBO components with or without scaling X.
#
# #' @importFrom Rcpp evalCpp
# #' @importFrom RcppParallel RcppParallelLibs
# #' @useDynLib mr.mash.alpha
# #'
# inner_loop_general_Rcpp <- function(n, XtXmu1, mu1, V, Vinv, w0, S0, precomp_quants,
#                                     standardize, compute_ELBO, update_V, update_order,
#                                     eps, nthreads){
# 
#   out <- inner_loop_general_rcpp(n, XtXmu1, mu1, V, Vinv, w0, S0, precomp_quants,
#                                  standardize, compute_ELBO, update_V, update_order,
#                                  eps, nthreads)
# 
#   ###Return output
#   if(compute_ELBO && update_V){
#     return(list(mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_tr_wERSS=out$var_part_tr_wERSS,
#                 neg_KL=out$neg_KL, var_part_ERSS=out$var_part_ERSS))
#   } else if(compute_ELBO && !update_V){
#     return(list(mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_tr_wERSS=out$var_part_tr_wERSS,
#                 neg_KL=out$neg_KL))
#   } else if(!compute_ELBO && update_V) {
#     return(list(mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_ERSS=out$var_part_ERSS))
#   } else {
#     return(list(mu1=out$mu1, S1=out$S1, w1=out$w1))
#   }
# }

###Wrapper of the inner loop with R or Rcpp
inner_loop_general_rss <- function(n, XtY, XtXmu1, mu1, V, Vinv, w0, S0, precomp_quants, 
                               standardize, compute_ELBO, update_V, version,
                               update_order, eps, nthreads){
  if(version=="R"){
    out <- inner_loop_general_rss_R(n, XtY, XtXmu1, mu1, V, Vinv, w0, S0, precomp_quants, 
                                standardize, compute_ELBO, update_V, update_order, eps)
  } # else if(version=="Rcpp"){
  #   update_order <- as.integer(update_order-1)
  #   out <- inner_loop_general_Rcpp(n, XtXmu1, mu1, V, Vinv, w0, simplify2array_custom(S0), precomp_quants, 
  #                                  standardize, compute_ELBO, update_V, update_order, eps, nthreads)
  # }
  # 
  return(out)
}


###Perform one iteration of the outer loop with or without scaling X
mr_mash_update_general_rss <- function(n, XtX, XtY, YtY, mu1_t, V, Vinv, ldetV, w0, S0,
                                      precomp_quants, compute_ELBO, standardize, 
                                      update_V, version, update_order, eps, 
                                      nthreads){
  
  

  ##Compute ??
  XtXmu1 <- XtX%*%mu1_t
  
  ##Update variational parameters, expected residuals, and ELBO components
  updates <- inner_loop_general_rss(n=n, XtY=XtY, XtXmu1=XtXmu1, mu1=mu1_t, V=V, Vinv=Vinv, w0=w0, S0=S0, 
                                precomp_quants=precomp_quants, standardize=standardize,
                                compute_ELBO=compute_ELBO, update_V=update_V, version=version,
                                update_order=update_order, eps=eps, nthreads=nthreads)   
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  
  out <- list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t)
  
  if(compute_ELBO || update_V){
    RbartRbar <- YtY - crossprod(mu1_t, XtY) - crossprod(XtY, mu1_t) + crossprod(mu1_t, XtX)%*%mu1_t
  }
  
  if(compute_ELBO && update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    out$ELBO <- compute_ELBO_rss_fun(n=n, RbartRbar=RbartRbar, Vinv=Vinv, ldetV=ldetV, 
                                     var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
    out$var_part_ERSS <- updates$var_part_ERSS
    out$RbartRbar
    
  } else if(compute_ELBO && !update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    out$ELBO <- compute_ELBO_rss_fun(n=n, RbartRbar=RbartRbar, Vinv=Vinv, ldetV=ldetV, 
                                     var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL)
    
  } else if(!compute_ELBO && update_V){
    out$var_part_ERSS <- updates$var_part_ERSS
    out$RbartRbar
  }
  
  return(out)
}


###Update V
update_V_rss_fun <- function(n, RbartRbar, var_part_ERSS){

  ERSS <- RbartRbar + var_part_ERSS
  V <- ERSS/n
  
  return(V)
}


###Update mixture weights
update_weights_em <- function(x){
  w <- colSums(x)
  w <- w/sum(w)
  return(w)
}
