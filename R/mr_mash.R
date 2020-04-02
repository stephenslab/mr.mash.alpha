#' @title  Multiple Regression with Multivariate Adaptive Shrinkage.
#' 
#' @description Performs multivariate multiple regression with
#'   mixture-of-normals prior.
#' 
#' @param Y n x r matrix of responses.
#' 
#' @param X n x p matrix of covariates.
#' 
#' @param V r x r residual covariance matrix.
#' 
#' @param S0 List of length K containing the desired r x r prior
#'   covariance matrices on the regression coefficients.
#' 
#' @param w0 K-vector with prior mixture weights, each associated with
#'   the respective covariance matrix in \code{S0}.
#' 
#' @param mu1_init p x r matrix of initial estimates of the posterior
#'   mean regression coefficients.
#' 
#' @param tol Convergence tolerance.
#' 
#' @param max_iter Maximum number of iterations for the optimization
#'   algorithm.
#' 
#' @param update_w0 If \code{TRUE}, prior weights are updated.
#' 
#' @param update_w0_method Method to update prior weights.
#' 
#' @param update_V if \code{TRUE}, residual covariance is updated.
#' 
#' @param compute_ELBO If \code{TRUE}, ELBO is computed.
#' 
#' @param standardize If \code{TRUE}, X is "standardized" using the
#'   sample means and sample standard deviations. Standardizing X
#'   allows a faster implementation, but the prior has a different
#'   interpretation. Coefficients and covariances are returned on the
#'   original scale.
#' 
#' @param version Whether to use R or C++ code to perform the
#'   coordinate ascent updates.
#' 
#' @param verbose If \code{TRUE}, some information about the
#'   algorithm's process is printed at each iteration.
#' 
#' @param e A small number to add to the diagonal elements of the
#'   prior matrices to improve numerical stability of the updates.
#'   
#' @param ca_update_order The order with which coordinated are updated.
#'   So far, only "consecutive" is supported.
#' 
#' @return A mr.mash fit, stored as a list with some or all of the
#' following elements:
#' 
#' \item{mu1}{p x r matrix of posterior means for the regression
#'   coeffcients.}
#' 
#' \item{S1}{r x r x p array of posterior covariances for the
#'   regression coeffcients.}
#' 
#' \item{w1}{p x K matrix of posterior assignment probabilities to the
#'   mixture components.}
#'   
#' \item{V}{r x r residual covariance matrix}
#' 
#' \item{w0}{K-vector with (updated, if \code{update_w0=TRUE}) prior mixture weights, each associated with
#'   the respective covariance matrix in \code{S0}}.
#'   
#' \item{S0}{r x r x K array of prior covariance matrices
#'   on the regression coefficients}.
#' 
#' \item{intercept}{r-vector containing posterior mean estimate of the
#'   intercept.}
#' 
#' \item{fitted}{n x r matrix of fitted values.}
#' 
#' \item{ELBO}{Evidence Lower Bound (ELBO) at last iteration.}
#' 
#' \item{progress}{A data frame including information regarding
#'   convergence criteria at each iteration.}
#'  
#' 
#' @examples 
#' ###Set seed
#' set.seed(123)
#'
#' ###Simulate X and Y
#' ##Set parameters
#' n  <- 100
#' p <- 10
#'
#' ##Compute residual covariance
#' V  <- rbind(c(1.0,0.2),
#'             c(0.2,0.4))
#'
#' ##Set true effects
#' B  <- matrix(c(-2, -2,
#'                5, 5,
#'                rep(0, (p-2)*2)), byrow=TRUE, ncol=2)
#'
#' ##Simulate X
#' X <- matrix(rnorm(n*p), nrow=n, ncol=p)
#' X <- scale(X, center=TRUE, scale=FALSE)
#'
#' ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix
#' ##and V is the residual covariance  
#' Y <- mr.mash.alpha:::sim_mvr(X, B, V)
#'
#' ###Specify the mixture weights and covariance matrices for the
#' ### mixture-of-normals prior.
#' grid <- seq(1, 5)
#' S0mix <- mr.mash.alpha:::compute_cov_canonical(ncol(Y), singletons=TRUE,
#'            hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid, zeromat=TRUE)
#' w0    <- rep(1/(length(S0mix)), length(S0mix))
#' 
#' ###Split the data in training and test sets
#' Ytrain <- Y[-c(1:10), ]
#' Xtrain <- X[-c(1:10), ]
#' Ytest <- Y[c(1:10), ]
#' Xtest <- X[c(1:10), ]
#'
#' ###Estimate residual covariance
#' V_est <- cov(Ytrain)
#'
#' ###Fit mr.mash
#' fit <- mr.mash(Xtrain, Ytrain, S0mix, w0, V_est, tol=1e-8, update_w0=TRUE,
#'                update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE,
#'                verbose=TRUE, update_V=TRUE, version="R", e=1e-8)
#'
#' # Compare the "fitted" values of Y against the true Y in the training set.
#' plot(fit$fitted,Ytrain,pch = 20,col = "darkblue",xlab = "true",
#'      ylab = "fitted")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted")
#'
#' # Predict the multivariate outcomes in the test set using the fitted model.
#' Ytest_est <- predict(fit,Xtest)
#' plot(Ytest_est,Ytest,pch = 20,col = "darkblue",xlab = "true",
#'      ylab = "predicted")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted")
#' 
#' @export
#' 
mr.mash <- function(X, Y, S0, w0=rep(1/(length(S0)), length(S0)), V=cov(Y), 
                    mu1_init=matrix(0, nrow=ncol(X), ncol=ncol(Y)), tol=1e-4,
                    max_iter=5000, update_w0=TRUE, update_w0_method=c("EM", "mixsqp"), 
                    compute_ELBO=TRUE, standardize=TRUE, verbose=TRUE,
                    update_V=FALSE, version=c("Rcpp", "R"), e=1e-8,
                    ca_update_order=c("consecutive", "decreasing_logBF")) {

  tic <- Sys.time()
  cat("Processing the inputs... ")

  # CHECK AND PROCESS INPUTS
  # ------------------------
  ###Select method to update the weights (if not specified by user, EM
  ###will be used)
  update_w0_method <- match.arg(update_w0_method)
  
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
  if(sum(w0)!=1)
    stop("Elements of w0 must sum to 1.")
  if(length(S0)!=length(w0))
    stop("S0 and w0 must have the same length.")
  if(update_w0_method=="mixsqp" && !compute_ELBO)
    stop("ELBO needs to be computed with update_w0_method=\"mixsqp\".")

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
  
  ###Initilize mu1, S1, w1, delta_mu1, ELBO, iterator, and progress
  mu1_t <- mu1_init 
  delta_mu1 <- matrix(Inf, nrow=p, ncol=r)
  ELBO <- -Inf
  progress <- as.data.frame(matrix(NA, nrow=max_iter, ncol=2))
  colnames(progress) <- c("iter", "mu1_max.diff")
  if(compute_ELBO){
    progress$ELBO_diff <- NA
    progress$ELBO <- NA
  }
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
    update_order <- compute_rank_variables_BFmix(X, Y, V, Vinv, w0, S0, comps, standardize, version)
  }
  
  cat("Done!\n")

  # PERFORM ONE UPDATE
  # ------------------
  ###First iteration
  t <- 0
  cat("Fitting the optimization algorithm... ")
  if(verbose){
    cat("\n")
    cat(" iter    mu1_max.diff")
    if(compute_ELBO)
      cat("     ELBO_diff               ELBO\n")
    else 
      cat("\n")
  }
  
  ##Save current estimates.
  mu1_tminus1 <- mu1_t   
  
  ##Update iterator
  t <- t+1
  
  ##Set last value of ELBO as ELBO0
  ELBO0 <- ELBO
  
  ###Update variational parameters
  ups <- mr_mash_update_general(X=X, Y=Y, mu1_t=mu1_t, V=V, Vinv=Vinv,
                                ldetV=ldetV, w0=w0, S0=S0,
                                precomp_quants=comps,
                                compute_ELBO=compute_ELBO,
                                standardize=standardize, 
                                update_V=update_V, version=version,
                                update_order=update_order)
  mu1_t <- ups$mu1_t
  S1_t  <- ups$S1_t
  w1_t  <- ups$w1_t
  if(compute_ELBO)
    ELBO <- ups$ELBO
  if(update_V)
    var_part_ERSS <- ups$var_part_ERSS
  
  ##Update progress data.frame 
  progress[t, c(1, 2)] <- c(t, max(delta_mu1))
  if(compute_ELBO)
    progress[t, c(3, 4)] <- c(ELBO - ELBO0, ELBO)
  
  if(verbose){
    ##Print out useful info
    cat(sprintf("%4d      %9.2e", t, max(delta_mu1)))
    if(compute_ELBO)
      cat(sprintf("      %9.2e      %0.20e\n", ELBO - ELBO0, ELBO))
    else
      cat("\n")
  }
  
  # MAIN LOOP
  # ---------
  ###Repeat the following until convergence, or until maximum number
  ###of iterations is reached.
  while(any(delta_mu1>tol)){
      
    ##Save current estimates.
    mu1_tminus1 <- mu1_t   
    
    ##Update iterator
    t <- t+1

    ##Exit loop if maximum number of iterations is reached
    if(t>max_iter){
      warning("Max number of iterations reached. Try increasing max_iter.")
      break
    }
    
    ##Set last value of ELBO as ELBO0
    ELBO0 <- ELBO

    # M-STEP
    # ------
    ##Update V if requested
    if(update_V){
      V     <- update_V_fun(Y, X, mu1_t, var_part_ERSS)
      
      #Recompute precomputed quantities after updating V
      comps <- precompute_quants(X, V, S0, standardize, version)
      if(!standardize)
        comps$xtx <- xtx
      if(compute_ELBO || !standardize)
        Vinv <- chol2inv(comps$V_chol)
      if(compute_ELBO)
        ldetV <- chol2ldet(comps$V_chol)
    }
    
    ##Update w0 if requested
    if(update_w0){
      if(update_w0_method=="EM")
        w0 <- update_weights_em(w1_t)
      else if(update_w0_method=="mixsqp"){
        w0em <- update_weights_em(w1_t)
        w0   <- update_weights_mixsqp(X=X, Y=Y, mu1=mu1_t, V=V, Vinv=Vinv,
                                      ldetV=ldetV, w0em=w0em, S0=S0,
                                      precomp_quants=comps,
                                      standardize=standardize,
                                      version=version, update_order=update_order)$w0
      }
    }

    # E-STEP
    # ------
    ###Update variational parameters
    ups <- mr_mash_update_general(X=X, Y=Y, mu1_t=mu1_t, V=V, Vinv=Vinv,
                                  ldetV=ldetV, w0=w0, S0=S0, 
                                  precomp_quants=comps,
                                  compute_ELBO=compute_ELBO,
                                  standardize=standardize,
                                  update_V=update_V, version=version, 
                                  update_order=update_order)
    mu1_t <- ups$mu1_t
    S1_t  <- ups$S1_t
    w1_t  <- ups$w1_t
    if(compute_ELBO)
      ELBO <- ups$ELBO
    if(update_V)
      var_part_ERSS <- ups$var_part_ERSS
    
    ##Compute distance in mu1 between two successive iterations
    delta_mu1 <- abs(mu1_t - mu1_tminus1)
    
    ##Update progress data.frame 
    progress[t, c(1, 2)] <- c(t, max(delta_mu1))
    if(compute_ELBO)
      progress[t, c(3, 4)] <- c(ELBO - ELBO0, ELBO)
    
    if(verbose){
      ##Print out useful info
      cat(sprintf("%4d      %9.2e", t, max(delta_mu1)))
      if(compute_ELBO)
        cat(sprintf("      %9.2e      %0.20e\n", ELBO - ELBO0, ELBO))
      else
        cat("\n")
    }
  }
  
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
  dimnames(S1_t)[[1]] <- colnames(Y)
  dimnames(S1_t)[[2]] <- colnames(Y)
  dimnames(S1_t)[[3]] <- colnames(X)
  rownames(w1_t) <- colnames(X)
  colnames(w1_t) <- names(S0)
  rownames(V) <- colnames(Y)
  colnames(V) <- colnames(Y)
  rownames(fitted_vals) <- rownames(Y)
  colnames(fitted_vals) <- colnames(Y)
  
  ###Remove unused rows of progress
  progress <- progress[rowSums(is.na(progress)) != ncol(progress), ]
  
  ###Return the posterior assignment probabilities (w1), the
  ###posterior mean of the coefficients (mu1), and the posterior
  ###covariance of the coefficients (S1), the residual covariance (V),
  ###the prior weights (w0), the intercept (intercept), the fitted values (fitted), 
  ###and the progress data frame (progress), and the prior covariance (S0) and,
  ###if computed, the Evidence Lower Bound (ELBO).
  out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, V=V, w0=w0, S0=simplify2array_custom(S0), 
              intercept=intercept, fitted=fitted_vals, progress=progress)
  if(compute_ELBO)
    ###Append ELBO to the output
    out$ELBO <- ELBO

  class(out) <- c("mr.mash", "list")
  
  cat("Done!\n")
  toc <- Sys.time()
  cat("mr.mash successfully executed in", difftime(toc,tic, units="mins"),
      "minutes!\n")
  
  return(out)
}
