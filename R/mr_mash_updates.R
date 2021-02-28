###Update variational parameters, expected residuals, and ELBO components with or without scaling X
inner_loop_general_R <- function(X, Rbar, mu1, V, Vinv, w0, S0, ###note: V is only needed when not scaling X
                               precomp_quants, standardize, compute_ELBO, update_V,
                               update_order, eps){
  ###Create variables to store quantities
  n <- nrow(Rbar)
  r <- ncol(Rbar)
  p <- ncol(X)
  K <- length(S0)
  S1    <- array(0, c(r, r, p))
  w1    <- matrix(0, nrow=p, ncol=K)
  var_part_tr_wERSS <- 0
  neg_KL <- 0
  var_part_ERSS <- matrix(0, nrow=r, ncol=r)
  
  ##Loop through the variables
  for(j in update_order){
    
    #Remove j-th effect from expected residuals 
    Rbar_j <- Rbar + outer(X[, j], mu1[j, ])
    
    #Run Bayesian SLR
    if(standardize){
      bfit <- bayes_mvr_mix_standardized_X(X[, j], Rbar_j, w0, S0, precomp_quants$S, precomp_quants$S1, 
                                     precomp_quants$SplusS0_chol, precomp_quants$S_chol, eps)      
    } else {
      bfit <- bayes_mvr_mix_centered_X(X[, j], Rbar_j, V, w0, S0, precomp_quants$xtx[j], Vinv, 
                                          precomp_quants$V_chol, precomp_quants$d, 
                                          precomp_quants$QtimesV_chol, eps)
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
      ELBO_parts <- compute_ELBO_terms(var_part_tr_wERSS, neg_KL, X[, j], Rbar_j, bfit, xtx, Vinv)
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
    Rbar <- Rbar_j - outer(X[, j], mu1[j, ])
  }
  
  ###Return output
  if(compute_ELBO && update_V){
    return(list(Rbar=Rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL, var_part_ERSS=var_part_ERSS))
  } else if(compute_ELBO && !update_V){
    return(list(Rbar=Rbar, mu1=mu1, S1=S1, w1=w1, var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
  } else if(!compute_ELBO && update_V) {
    return(list(Rbar=Rbar, mu1=mu1, S1=S1, w1=w1, var_part_ERSS=var_part_ERSS))
  } else { 
    return(list(Rbar=Rbar, mu1=mu1, S1=S1, w1=w1))
  }
}

### Wrapper for the Rcpp function to update variational parameters,
### expected residuals, and ELBO components with or without scaling X.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib mr.mash.alpha
#' 
inner_loop_general_Rcpp <- function(X, Rbar, mu1, V, Vinv, w0, S0, precomp_quants, 
                                    standardize, compute_ELBO, update_V, update_order, 
                                    eps, nthreads){
  
  out <- inner_loop_general_rcpp(X, Rbar, mu1, V, Vinv, w0, S0, precomp_quants, 
                                 standardize, compute_ELBO, update_V, update_order,
                                 eps, nthreads)
  
  ###Return output
  if(compute_ELBO && update_V){
    return(list(Rbar=out$Rbar, mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_tr_wERSS=out$var_part_tr_wERSS, 
                neg_KL=out$neg_KL, var_part_ERSS=out$var_part_ERSS))
  } else if(compute_ELBO && !update_V){
    return(list(Rbar=out$Rbar, mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_tr_wERSS=out$var_part_tr_wERSS, 
                neg_KL=out$neg_KL))
  } else if(!compute_ELBO && update_V) {
    return(list(Rbar=out$Rbar, mu1=out$mu1, S1=out$S1, w1=out$w1, var_part_ERSS=out$var_part_ERSS))
  } else { 
    return(list(Rbar=out$Rbar, mu1=out$mu1, S1=out$S1, w1=out$w1))
  }
}

###Wrapper of the inner loop with R or Rcpp
inner_loop_general <- function(X, Rbar, mu1, V, Vinv, w0, S0, precomp_quants, 
                               standardize, compute_ELBO, update_V, version,
                               update_order, eps, nthreads){
  if(version=="R"){
    out <- inner_loop_general_R(X, Rbar, mu1, V, Vinv, w0, S0, precomp_quants, 
                                standardize, compute_ELBO, update_V, update_order, eps)
  } else if(version=="Rcpp"){
    update_order <- as.integer(update_order-1)
    out <- inner_loop_general_Rcpp(X, Rbar, mu1, V, Vinv, w0, simplify2array_custom(S0), precomp_quants, 
                                   standardize, compute_ELBO, update_V, update_order, eps, nthreads)
  }
  
  return(out)
}


###Perform one iteration of the outer loop with or without scaling X
mr_mash_update_general <- function(X, Y, mu1_t, mu, V, Vinv, ldetV, w0, S0,
                                   precomp_quants, compute_ELBO, standardize, 
                                   update_V, version, update_order, eps, 
                                   nthreads, Y_miss_patterns){
  
  
  if(!is.null(Y_miss_patterns)){
    ##Impute missing Ys
    outY <- impute_missing_Y(Y=Y, mu=mu, Vinv=Vinv, miss=Y_miss_patterns$miss, non_miss=Y_miss_patterns$non_miss, 
                             version=version)
    Y <- outY$Y
    Y_cov <- outY$Y_cov
    sum_neg_ent_Y_miss <- outY$sum_neg_ent_Y_miss
    
    # Update the intercept
    muy <- colMeans(Y)
    
    ##Compute expected residuals
    Rbar <- scale_fast2(Y, scale=FALSE)$M - X%*%mu1_t
    
  } else {
    Y_cov <- NULL
    sum_neg_ent_Y_miss <- 0
    
    ##Compute expected residuals
    Rbar <- Y - mu
  }
  
  ##Update variational parameters, expected residuals, and ELBO components
  updates <- inner_loop_general(X=X, Rbar=Rbar, mu1=mu1_t, V=V, Vinv=Vinv, w0=w0, S0=S0, 
                                precomp_quants=precomp_quants, standardize=standardize,
                                compute_ELBO=compute_ELBO, update_V=update_V, version=version,
                                update_order=update_order, eps=eps, nthreads=nthreads)   
  mu1_t   <- updates$mu1
  S1_t    <- updates$S1
  w1_t    <- updates$w1
  Rbar    <- updates$Rbar
  
  out <- list(mu1_t=mu1_t, S1_t=S1_t, w1_t=w1_t)
  
  if(!is.null(Y_miss_patterns)){
    out$Y <- Y 
    out$muy <- muy
    out$Y_cov <- Y_cov
  }
  
  if(compute_ELBO && update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    out$ELBO <- compute_ELBO_fun(Rbar=Rbar, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS, 
                                 neg_KL=neg_KL, Y_cov=Y_cov, sum_neg_ent_Y_miss=sum_neg_ent_Y_miss)
    
    out$var_part_ERSS <- updates$var_part_ERSS
    
  } else if(compute_ELBO && !update_V){
    ##Compute ELBO
    var_part_tr_wERSS <- updates$var_part_tr_wERSS
    neg_KL <- updates$neg_KL
    out$ELBO <- compute_ELBO_fun(Rbar=Rbar, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=var_part_tr_wERSS,
                                 neg_KL=neg_KL, Y_cov=Y_cov, sum_neg_ent_Y_miss=sum_neg_ent_Y_miss)
    
  } else if(!compute_ELBO && update_V){
    out$var_part_ERSS <- updates$var_part_ERSS
  }
  
  return(out)
}


###Update V
update_V_fun <- function(Y, mu, var_part_ERSS, Y_cov){
  n <- nrow(Y)
  
  Rbar <- Y - mu
  ERSS <- crossprod(Rbar) + var_part_ERSS + Y_cov
  V <- ERSS/n
  
  return(V)
}


###Update mixture weights
update_weights_em <- function(x){
  w <- colSums(x)
  w <- w/sum(w)
  return(w)
}


###Impute/update missing Y
impute_missing_Y_R <- function(Y, mu, Vinv, miss, non_miss){
  n <- nrow(Y)
  r <- ncol(Y)
  
  Y_cov <- matrix(0, r, r)
  sum_neg_ent_Y_miss <- 0
  
  for (i in 1:n){
    non_miss_i <- non_miss[[i]]
    miss_i <- miss[[i]]
    Vinv_mo <- Vinv[miss_i, non_miss_i, drop=FALSE]
    Vinv_mm <- Vinv[miss_i, miss_i, drop=FALSE]
    if(any(miss_i)){
      # Compute variance
      Y_cov_i <- matrix(0, r, r)
      Vinv_mm_chol <- chol(Vinv_mm)
      Y_cov_mm <- chol2inv(Vinv_mm_chol)
      Y_cov_i[miss_i, miss_i] <- Y_cov_mm
      
      Y_cov <- Y_cov + Y_cov_i
      
      # Compute mean
      Y[i, miss_i] <- mu[i, miss_i] - Y_cov_mm %*% Vinv_mo %*% (Y[i, non_miss_i] - mu[i, non_miss_i])
      
      # Compute sum of the negative entropy of Y missing
      #sum_neg_ent_Y_miss <- sum_neg_ent_Y_miss + (0.5 * as.numeric(determinant(1/(2*pi*exp(1))*Vinv_mm, logarithm = TRUE)$modulus))
      sum_neg_ent_Y_miss <- sum_neg_ent_Y_miss + (0.5 * (ncol(Vinv_mm_chol)*log((1/(2*pi*exp(1)))) + chol2ldet(Vinv_mm_chol))) # log(det(kA)) = r*log(k) + log(det(A)) where is the size of the matrix
    }
  }
  
  return(list(Y=Y, Y_cov=Y_cov, sum_neg_ent_Y_miss=sum_neg_ent_Y_miss))
}

###Wrapper of impute/update missing Y with R or Rcpp
impute_missing_Y <- function(Y, mu, Vinv, miss, non_miss, version){
  if(version=="R"){
    out <- impute_missing_Y_R(Y, mu, Vinv, miss, non_miss)
  } else if(version=="Rcpp"){
    out <- impute_missing_Y_rcpp(Y, mu, Vinv, simplify2array_custom(miss), simplify2array_custom(non_miss))
  }
  
  return(out)
}