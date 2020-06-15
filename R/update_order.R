###Compute logbf from Bayesian multivariate simple regression with mixture prior
compute_logbf_R <- function(X, Y, V, Vinv, w0, S0, precomp_quants, standardize, eps){
  p <- ncol(X)
  logbf <- rep(0, p)
  
  for(j in 1:p){
    #Run Bayesian SLR
    if(standardize){
      bfit <- bayes_mvr_mix_standardized_X(X[, j], Y, w0, S0, precomp_quants$S, precomp_quants$S1, 
                                           precomp_quants$SplusS0_chol, precomp_quants$S_chol, eps)      
    } else {
      bfit <- bayes_mvr_mix_centered_X(X[, j], Y, V, w0, S0, precomp_quants$xtx[j], Vinv, 
                                       precomp_quants$V_chol, precomp_quants$d, 
                                       precomp_quants$QtimesV_chol, eps)
    }
    
    logbf[j] <- bfit$logbf
  }
  
  return(logbf)
}

###Compute rank of logbf from Bayesian multivariate simple regression with mixture prior
compute_rank_variables_BFmix <- function(X, Y, V, Vinv, w0, S0, precomp_quants, standardize, version, decreasing, eps){
  if(version=="R"){
    logbfs <- compute_logbf_R(X, Y, V, Vinv, w0, S0, precomp_quants, standardize, eps)
  } else if(version=="Rcpp"){
    logbfs <- compute_logbf_rcpp(X, Y, V, Vinv, w0, simplify2array_custom(S0), precomp_quants, standardize, eps)
    logbfs <- drop(logbfs)
  }
  
  if(decreasing){
    ##Negative sign is needed because rank() by default ranks from smallest to largest
    ##while we want from largest to smallest
    rank_variables_BFmix <- rank(-logbfs, ties.method="first", na.last="keep")
  } else {
    rank_variables_BFmix <- rank(logbfs, ties.method="first", na.last="keep")
  }
  
  return(rank_variables_BFmix)
}

