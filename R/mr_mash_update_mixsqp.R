###Update mixture weights with mixsqp
#' @importFrom mixsqp mixsqp
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib mr.mash.alpha
compute_mixsqp_update <- function (X, Y, V, S0, mu1, Vinv, precomp_quants,
                                   standardize, version, update_order) {
  
  # Compute the p x k matrix of log-likelihoods conditional on each
  # prior mixture component.
  Rbar <- Y - X%*%mu1
  L <- compute_mixsqp_update_loop_general(X, Rbar, V, S0, mu1, Vinv, precomp_quants, standardize, version, update_order)
  
  # Run mixsqp
  out <- mixsqp(L,log = TRUE,control = list(verbose = FALSE))

  # Return the updated mixture weights ("w0") and the number of
  # mix-SQP iterations performed.
  return(list(w0=out$x, numiter=nrow(out$progress)))
}

###Wrapper of the loop to compute mixsqp update of the weights
compute_mixsqp_update_loop_R <- function(X, Rbar, V, S0, mu1, Vinv, precomp_quants, standardize, update_order){
  # Get the number of predictors (p), the number of mixture
  # components in the prior (k), and the number of samples (n).
  p <- ncol(X)
  K <- length(S0)
  n <- nrow(Rbar)
  
  L <- matrix(0, p, K)
  
  for(j in update_order){
    #Remove j-th effect from expected residuals 
    Rbar_j <- Rbar + outer(X[, j], mu1[j, ])
    
    if(standardize){
      # Compute the least-squares estimate.
      xtx <- n-1
      b <- drop(X[, j] %*% Rbar_j)/xtx
    } else {
      # Compute the least-squares estimate and covariance.
      b <- drop(X[, j] %*% Rbar_j)/precomp_quants$xtx[j]
      S <- V/precomp_quants$xtx[j]
      
      # Compute quantities needed for bayes_mvr_ridge_centered_X()
      S_chol <- precomp_quants$V_chol/sqrt(precomp_quants$xtx[j])
    }
    
    for(k in 1:K){
      if(standardize){
        L[j, k] <- bayes_mvr_ridge_standardized_X(b, S0[[k]], precomp_quants$S, 
                                            precomp_quants$S1[[k]], precomp_quants$SplusS0_chol[[k]], 
                                            precomp_quants$S_chol)$logbf        
      } else {
        L[j, k] <- bayes_mvr_ridge_centered_X(V, b, S, S0[[k]], precomp_quants$xtx[j], Vinv,
                                              precomp_quants$V_chol, S_chol, precomp_quants$d[[k]], 
                                              precomp_quants$QtimesV_chol[[k]])$logbf
      }
    }
  }
  
  return(L)
}

###Wrapper of the loop to compute mixsqp update of the weights with R or Rcpp
compute_mixsqp_update_loop_general <- function(X, Rbar, V, S0, mu1, Vinv, precomp_quants, standardize, version,
                                               update_order){
  if(version=="R"){
    out <- compute_mixsqp_update_loop_R(X, Rbar, V, S0, mu1, Vinv, precomp_quants, standardize, update_order)
  } else if(version=="Rcpp"){
    update_order <- update_order-1
    out <- compute_mixsqp_update_loop_rcpp(X, Rbar, V, simplify2array_custom(S0), mu1, Vinv, precomp_quants, standardize,
                                           update_order)
  }
  
  return(out)
}

# Perform backtracking line search to identify a step size for the
# mixture weights update that increases the ELBO.
backtracking_line_search <- function (X, Y, V, Vinv, ldetV, S0, mu1, w0old, w0mixsqp,
                                      precomp_quants, standardize, compute_ELBO, update_V, 
                                      version, update_order, stepsize.reduce, stepsize.min) {
  
  # Compute the objective (ELBO) at the current iterate.
  ##Update variational parameters, expected residuals, and ELBO components
  Rbar <- Y - X%*%mu1
  
  updates <- inner_loop_general(X=X, Rbar=Rbar, mu1=mu1, V=V, Vinv=Vinv, w0=w0old, S0=S0, 
                                precomp_quants=precomp_quants, standardize=standardize, 
                                compute_ELBO=compute_ELBO, update_V=update_V, version=version,
                                update_order=update_order)

  f <- compute_ELBO_fun(Rbar=Rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=updates$var_part_tr_wERSS, neg_KL=updates$neg_KL)
    
  
  # Perform backtracking line search to identify a step size that
  # increases the ELBO.
  a <- 1
  while (TRUE) {
    w0new <- a*w0mixsqp + (1 - a)*w0old
    updates <- inner_loop_general(X=X, Rbar=Rbar, mu1=mu1, V=V, Vinv=Vinv, w0=w0new, S0=S0, 
                                  precomp_quants=precomp_quants, standardize=standardize,
                                  compute_ELBO=compute_ELBO, update_V=update_V, version=version,
                                  update_order=update_order) 

    fnew <- compute_ELBO_fun(Rbar=Rbar, V=V, Vinv=Vinv, ldetV=ldetV, var_part_tr_wERSS=updates$var_part_tr_wERSS, neg_KL=updates$neg_KL)
    
    # Check whether the new candidate increases the ELBO.
    if (fnew > f)
      break
    
    # If we cannot decrease the step size further, terminate the
    # backtracking line search, and set the step size to be the
    # minimum step size.
    else if (a * stepsize.reduce < stepsize.min) {
      
      # We need to terminate backtracking line search because we have
      # arrived at the smallest allowable step size.
      a     <- 0
      w0new <- w0old
      break
    }
    
    # The new candidate does not increase the ELBO, so we need to try
    # again with a smaller step size.
    a = a * stepsize.reduce
  }
  
  # Return the updated mixture weights ("w0") and the step size
  # determined by the backtracking line search ("a").
  return(list(w0 = w0new,a = a))
}

# Update the mixture weights with mix-SQP, following by backtracking
# line search to ensure that the ELBO does not decrease.
update_weights_mixsqp <- function (X, Y, mu1, V, Vinv, ldetV, w0old, S0,
                                   precomp_quants, standardize,
                                   compute_ELBO=TRUE, update_V=FALSE, version, update_order,
                                   stepsize.reduce = 0.5, stepsize.min = 1e-8) {
  
  # Compute the mix-SQP update for the mixture weights. Note that this
  # update is not guaranteed to increase the ELBO.
  out1 <- compute_mixsqp_update(X, Y, V, S0, mu1, Vinv, precomp_quants, standardize, version, update_order)
  
  # Perform backtracking line search to identify a step size that
  # increases the ELBO.
  out2 <- backtracking_line_search(X, Y, V, Vinv, ldetV, S0, mu1, w0old, out1$w0,
                                   precomp_quants, standardize, compute_ELBO, update_V, 
                                   version, update_order, stepsize.reduce, stepsize.min)
  
  # Return the updated mixture weights ("w0"), the number of mix-SQP
  # iterations performed ("numiter"), and the step size determined by the
  # backtracking line search ("a").
  return(list(w0      = out2$w0,
              numiter = out1$numiter,
              a       = out2$a))
}

