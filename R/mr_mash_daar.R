#' @export
#' 
mr.mash.daar <- function(X, Y, S0, w0=rep(1/(length(S0)), length(S0)), V=cov(Y), 
                        mu1_init=matrix(0, nrow=ncol(X), ncol=ncol(Y)), tol=1e-4,
                        max_iter=5000
                        compute_ELBO=TRUE, standardize=TRUE,
                        update_V=FALSE, version=c("Rcpp", "R"), e=1e-8,
                        ca_update_order=c("consecutive", "decreasing_logBF", "increasing_logBF")) {
  
  tic <- Sys.time()
  cat("Processing the inputs... ")
  
  # CHECK AND PROCESS INPUTS
  # ------------------------
  ###Select version of the inner loop (if not specified by user, Rcpp
  ###will be used)
  version <- match.arg(version)
  
  ###Select ordering of the coordinate ascent updates (if not specified by user,
  ###consecutive will be used
  ca_update_order <- match.arg(ca_update_order)
  
  ###Check that the inputs are in the correct format
  if(!is.matrix(Y))
    stop("Y must be a matrix.")
  if(!is.matrix(X))
    stop("X must be a matrix.")
  if(any(is.na(Y)))
    stop("Y must not contain missing values.")
  if(any(is.na(X)))
    stop("X must not contain missing values.")
  if(!is.matrix(V) || !isSymmetric(V))
    stop("V must be a symmetric matrix.")
  if(!is.list(S0))
    stop("S0 must be a list.")
  if(!is.vector(w0))
    stop("w0 must be a vector.")
  if(abs(sum(w0) - 1) > 1e-10)
    stop("Elements of w0 must sum to 1.")
  if(length(S0)!=length(w0))
    stop("S0 and w0 must have the same length.")
  if(!is.matrix(mu1_init))
    stop("mu1_init must be a matrix.")

  ###Obtain dimensions needed from inputs
  p <- ncol(X)
  n <- nrow(X)
  r <- ncol(Y)
  K <- length(S0)
  
  # PRE-PROCESSING STEPS
  # --------------------
  ###Add number to diagonal elements of the prior matrices (improves
  ###numerical stability)
  S0 <- lapply(S0, makePD, e=e)
  
  ###Center Y, and center (and, optionally, scale) X
  outY <- scale_fast2(Y, scale=FALSE)
  outX <- scale_fast2(X, scale=standardize)
  muy <- outY$means
  mux <- outX$means
  if (standardize)
    sx <- outX$sds
  Y <- outY$M
  rm(outY)
  X <- outX$M
  rm(outX)
  
  ###Precompute quantities
  comps <- precompute_quants(X, V, S0, standardize, version)
  if(!standardize){
    xtx <- colSums(X^2)
    comps$xtx <- xtx
  }
  
  if(compute_ELBO || !standardize)
    ###Compute inverse of V (needed for the ELBO and unstandardized X)
    Vinv <- chol2inv(comps$V_chol)
  else {
    if(version=="R")
      Vinv <- NULL
    else if(version=="Rcpp")
      Vinv <- matrix(0, nrow=r, ncol=r)
  }
  
  if(compute_ELBO)
    ###Compute log determinant of V (needed for the ELBO)
    ldetV <- chol2ldet(comps$V_chol)
  else
    ldetV <- NULL
  
  ###Set the ordering of the coordinate ascent updates
  if(ca_update_order=="consecutive"){
    update_order <- 1:p
  } else if(ca_update_order=="decreasing_logBF"){
    update_order <- compute_rank_variables_BFmix(X, Y, V, Vinv, w0, S0, comps, standardize, version, decreasing=TRUE)
  } else if(ca_update_order=="increasing_logBF"){
    update_order <- compute_rank_variables_BFmix(X, Y, V, Vinv, w0, S0, comps, standardize, version, decreasing=FALSE)
  }
  
  cat("Done!\n")
  
 
  # MAIN LOOP
  # ---------
  
  out <- suppressWarnings(
    daarem(mu1_init, mr_mash_update_general_mu1_daar, mr_mash_update_general_ELBO_daar,
           X, Y, V, Vinv, ldetV, w0, S0, comps, compute_ELBO, standardize, update_V, version, update_order,
           control = list(maxiter = max_iter, order = 10, tol = tol,
                          mon.tol = 0.01,kappa = 20,alpha = 1.1)))
  
  cat("Done!\n")
  cat("Processing the outputs... ")
  
  # POST-PROCESSING STEPS
  # --------------------
  ###Compute the "fitted" values.
  fitted_vals <- addtocols(X %*% mu1_t, muy)
  
  if(standardize){
    ###Rescale posterior means and covariance of coefficients. In the
    ###context of predicting Y, this rescaling is equivalent to
    ###rescaling each column j of a given matrix, Xnew, by sx[j].
    post_rescaled <- rescale_post_mean_covar_fast(mu1_t, S1_t, sx)
    mu1_t <- post_rescaled$mu1_orig
    S1_t <- post_rescaled$S1_orig
  }
  
  ###Compute posterior mean estimate of intercept. Note that when
  ###columns of X are standardized, the intercept should be computed
  ###with respect to the *rescaled* coefficients to recover the
  ###correct fitted values. This is why this is done after rescaling
  ###the coefficients above.
  intercept <- drop(muy - mux %*% mu1_t)
  
  ###Assign names to outputs dimensions
  rownames(mu1_t) <- colnames(X)
  colnames(mu1_t) <- colnames(Y)
  rownames(V) <- colnames(Y)
  colnames(V) <- colnames(Y)
  rownames(fitted_vals) <- rownames(Y)
  colnames(fitted_vals) <- colnames(Y)
  
  
  ###Return the posterior assignment probabilities (w1), the
  ###posterior mean of the coefficients (mu1), and the posterior
  ###covariance of the coefficients (S1), the residual covariance (V),
  ###the prior weights (w0), the intercept (intercept), the fitted values (fitted), 
  ###and the progress data frame (progress), the prior covariance (S0), convergence
  ### status and, if computed, the Evidence Lower Bound (ELBO).
  out <- list(mu1=mu1_t, S0=simplify2array_custom(S0), 
              intercept=intercept, fitted=fitted_vals)
  if(compute_ELBO)
    ###Append ELBO to the output
    out$ELBO <- ELBO
  
  class(out) <- c("mr.mash", "list")
  
  cat("Done!\n")
  toc <- Sys.time()
  cat("mr.mash successfully executed in", difftime(toc, tic, units="mins"),
      "minutes!\n")
  
  return(out)
}




###Perform one iteration of the outer loop with or without scaling X
mr_mash_update_general_mu1_daar <- function(mu1_t, X, Y, V, Vinv, ldetV, w0, S0,
                                            precomp_quants, compute_ELBO, standardize, 
                                            update_V, version, update_order){
  
  out <- mr_mash_update_general(X, Y, mu1_t, V, Vinv, ldetV, w0, S0,
                                precomp_quants, compute_ELBO, standardize, 
                                update_V, version, update_order)
  mu1 <- out$mu1_t
  
  return(mu1)
}

###Perform one iteration of the outer loop with or without scaling X
mr_mash_update_general_ELBO_daar <- function(mu1_t, X, Y, V, Vinv, ldetV, w0, S0,
                                             precomp_quants, compute_ELBO, standardize, 
                                             update_V, version, update_order){
  out <- mr_mash_update_general(X, Y, mu1_t, V, Vinv, ldetV, w0, S0,
                                precomp_quants, compute_ELBO, standardize, 
                                update_V, version, update_order)
  ELBO <- out$ELBO
  
  return(ELBO)
}