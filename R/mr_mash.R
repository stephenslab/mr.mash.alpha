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
#' @param mu_init p x r matrix of initial estimates of the posterior
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
#' \item{intercept}{r-vector with the estimated intercepts.}
#' 
#' \item{fitted}{n x r matrix of fitted values.}
#' 
#' \item{ELBO}{Evidence Lower Bound (ELBO) at last iteration.}
#' 
#' \item{progress}{A data frame including information regarding
#'   convergence criteria at each iteration.}
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
#' fit <- mr.mash(Xtrain, Ytrain, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE,
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
mr.mash <- function(X, Y, V=NULL, S0, w0, mu_init=NULL, tol=1e-8,
                    max_iter=1e5, update_w0=TRUE,
                    update_w0_method=c("EM", "mixsqp"), 
                    compute_ELBO=TRUE, standardize=TRUE, verbose=TRUE,
                    update_V=FALSE, version=c("R", "Rcpp"), e=1e-8) {

  tic <- Sys.time()
  cat("Processing the inputs... ")

  # CHECK AND PROCESS INPUTS
  # ------------------------
  ###Select method to update the weights (if not specified by user, EM
  ###will be used)
  update_w0_method <- match.arg(update_w0_method)
  
  ###Select version of the inner loop (if not specified by user, R
  ###will be used)
  version <- match.arg(version)
  
  ###Check that the inputs are in the correct format, and initialize
  ###any model parameters that are not specified.
  if(!is.matrix(Y))
    stop("Y must be a matrix.")
  if(!is.matrix(X))
    stop("X must be a matrix.")
  if(any(is.na(Y)))
    stop("Y must not contain missing values.")
  if(any(is.na(X)))
    stop("X must not contain missing values.")
  if(is.null(V))
    V <- cov(Y)
  else if(!is.matrix(V) || !isSymmetric(V))
    stop("V must be a symmetric matrix.")
  if(!is.list(S0))
    stop("S0 must be a list.")
  if(!is.vector(w0))
    stop("w0 must be a vector.")
  if(length(S0)!=length(w0))
    stop("S0 and w0 must have the same length.")
  if(update_w0_method=="mixsqp" && !compute_ELBO)
    stop("ELBO needs to be computed with update_w0_method=\"mixsqp\".")

  # If not specified, set the initial estimates of the posterior mean
  # regression coefficients.
  p <- ncol(X)
  n <- nrow(X)
  R <- ncol(Y)
  K <- length(S0)
  if(is.null(mu_init))
    mu_init <- matrix(0, nrow=p, ncol=R)

  # PRE-PROCESSING STEPS
  # --------------------
  ###Add number to diagonal elements of the prior matrices (improves
  ###numerical stability)
  S0 <- lapply(S0, makePD, e=e)
  
  ###Center Y, and center (and, optionally, scale) X
  Y   <- scale(Y, center=TRUE, scale=FALSE)
  X   <- scale(X, center=TRUE, scale=standardize)
  muy <- attr(Y,"scaled:center")
  mux <- attr(X,"scaled:center")
  if (standardize)
    sx <- attr(X,"scaled:scale")
  attr(X,"scaled:center") <- NULL
  attr(X,"scaled:scale")  <- NULL
  attr(Y,"scaled:center") <- NULL
  
  ###Initilize mu1, S1, w1, error, ELBO, iterator, and progress
  mu1_t    <- mu_init 
  err      <- matrix(Inf, nrow=p, ncol=R)
  progress <- data.frame() 
  if(compute_ELBO)
    ELBO <- -Inf
  
  ###Precompute quantities
  comps <- precompute_quants(n, X, V, S0, standardize, version)
  
  if(compute_ELBO || !standardize)
    ###Compute inverse of V (needed for the ELBO and unstandardized X)
    #Vinv <- chol2inv(comps$V_chol)
    Vinv <- backsolve(comps$V_chol, forwardsolve(t(comps$V_chol),
                                                 diag(nrow(comps$V_chol))))
  else {
    if(version=="R")
      Vinv <- NULL
    else if(version=="Rcpp")
      Vinv <- matrix(0, nrow=R, ncol=R)
  }
  
  if(compute_ELBO)
    ###Compute log determinant of V (needed for the ELBO)
    ldetV <- chol2ldet(comps$V_chol)
  else
    ldetV <- NULL
  
  cat("Done!\n")

  ###First iteration
  t <- 0
  if(verbose){
    cat("Fitting the optimization algorithm... \n")
    if(compute_ELBO)
      cat(" iter    mu1_max.diff     ELBO_diff               ELBO\n")
    else 
      cat(" iter    mu1_max.diff\n")
  } else 
    cat("Fitting the optimization algorithm... ")
  
  ##Save current estimates.
  mu1_tminus1 <- mu1_t   
  
  ##Update iterator
  t <- t+1
  
  ##Set last value of ELBO as ELBO0
  if(compute_ELBO)
    ELBO0 <- ELBO
  
  ###Update variational parameters
  ups <- mr_mash_update_general(X=X, Y=Y, mu1_t=mu1_t, V=V, Vinv=Vinv,
                                ldetV=ldetV, w0=w0, S0=S0,
                                precomp_quants=comps,
                                compute_ELBO=compute_ELBO,
                                standardize=standardize, 
                                update_V=update_V, version=version)
  mu1_t <- ups$mu1_t
  S1_t  <- ups$S1_t
  w1_t  <- ups$w1_t
  if(compute_ELBO)
    ELBO <- ups$ELBO
  if(update_V)
    var_part_ERSS <- ups$var_part_ERSS
  
  if(compute_ELBO){
    if(verbose)
      ##Print out useful info
      cat(sprintf("%4d      %9.2e      %9.2e      %0.20e\n",
                  t, max(err), ELBO - ELBO0, ELBO))

    ##Update progress data.frame 
    progress <- rbind(progress, c(t, max(err), ELBO - ELBO0, ELBO))
  } else {
    if(verbose)
      ##Print out useful info
      cat(sprintf("%4d      %9.2e\n", t, max(err)))
    
    ##Update progress data.frame 
    progress <- rbind(progress, c(t, max(err)))
  }
  
  # MAIN LOOP
  # ---------
  ###Repeat the following until convergence, or until maximum number
  ###of iterations is reached.
  while(any(err>tol)){
      
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
    if(compute_ELBO)
      ELBO0 <- ELBO
    
    ##Update V if requested
    if(update_V){
      V     <- update_V_fun(Y, X, mu1_t, var_part_ERSS)
      comps <- precompute_quants(n, X, V, S0, standardize, version)
      if(compute_ELBO || !standardize)
        Vinv <- backsolve(comps$V_chol,
                          forwardsolve(t(comps$V_chol),
                                       diag(nrow(comps$V_chol))))
      
      ###Compute log determinant of V (needed for the ELBO)
      if(compute_ELBO)
        ldetV <- chol2ldet(comps$V_chol)
    }
    
    ##Update w0 if requested
    if(update_w0){
      if(update_w0_method=="EM")
        w0 <- update_weights_em(w1_t)
      else if(update_w0_method=="mixsqp"){
        w0em <- update_weights_em(w1_t)
        w0   <- update_weights_mixsqp(X=X, Y=Y, mu1_t=mu1_t, V=V, Vinv=Vinv,
                                      ldetV=ldetV, w0em=w0em, S0=S0,
                                      precomp_quants=comps,
                                      standardize=standardize,
                                      version=version)$w0
      }
    }
    
    ###Variational parameters
    ups <- mr_mash_update_general(X=X, Y=Y, mu1_t=mu1_t, V=V, Vinv=Vinv,
                                  ldetV=ldetV, w0=w0, S0=S0, 
                                  precomp_quants=comps,
                                  compute_ELBO=compute_ELBO,
                                  standardize=standardize,
                                  update_V=update_V, version=version)
    mu1_t <- ups$mu1_t
    S1_t  <- ups$S1_t
    w1_t  <- ups$w1_t
    if(compute_ELBO)
      ELBO <- ups$ELBO
    if(update_V)
      var_part_ERSS <- ups$var_part_ERSS
    
    ##Compute distance in mu1 between two successive iterations
    err <- abs(mu1_t - mu1_tminus1)
    
    if(compute_ELBO){
      if(verbose)
        ##Print out useful info
        cat(sprintf("%4d      %9.2e      %9.2e      %0.20e\n",
                    t, max(err), ELBO - ELBO0, ELBO))

      ##Update progress data.frame 
      progress <- rbind(progress, c(t, max(err), ELBO - ELBO0, ELBO))
    } else {
      if(verbose)
        ##Print out useful info
        cat(sprintf("%4d      %9.2e\n", t, max(err)))

      ##Update progress data.frame 
      progress <- rbind(progress, c(t, max(err)))
    }
  }
  
  cat("Done!\n")
  cat("Processing the output... ")
  
  ###Compute fitted values
  fitted_vals <- X %*% mu1_t
  fitted_vals <- addtocols(fitted_vals,muy)

  if(standardize){
    ###Rescale posterior means and covariance of coefficients
    mu1_t <- mu1_t/sx
    for(j in 1:dim(S1_t)[3])
      S1_t[, , j] <- S1_t[, , j]/(sx[j]^2)
  }

  ###Compute intercept
  intercept <- muy - mux %*% mu1_t
  intercept <- drop(intercept)
  
  if(compute_ELBO && update_V){
    colnames(progress) <- c("iter", "mu1_max.diff", "ELBO_diff", "ELBO")
    
    ###Return the posterior assignment probabilities (w1), the
    ###posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the residual covariance
    ###(V), the intercept (intercept), the fitted values (fitted), the
    ###Evidence Lower Bound (ELBO), and the progress data frame
    ###(progress).
    out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, V=V, intercept=intercept,
                fitted=fitted_vals, ELBO=ELBO, progress=progress)
  } else if(compute_ELBO && !update_V){
    colnames(progress) <- c("iter", "mu1_max.diff", "ELBO_diff", "ELBO")
    
    ###Return the posterior assignment probabilities (w1), the
    ###posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the intercept (intercept),
    ###the fitted values (fitted), the Evidence Lower Bound (ELBO),
    ###and the progress data frame (progress).
    out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept,
                fitted=fitted_vals, ELBO=ELBO, progress=progress)
  } else if(!compute_ELBO && update_V){
    colnames(progress) <- c("iter", "mu1_max.diff")
    
    ###Return the posterior assignment probabilities (w1), the
    ###posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the residual covariance
    ###(V), the intercept (intercept), the fitted values (fitted), and
    ###the progress data frame (progress).
    out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, V=V, intercept=intercept,
                fitted=fitted_vals, progress=progress)
  } else {
    colnames(progress) <- c("iter", "mu1_max.diff")
    
    ###Return the posterior assignment probabilities (w1), the
    ###posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the intercept (intercept),
    ###the fitted values (fitted), and the progress data frame
    ###(progress).
    out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept,
                fitted=fitted_vals, progress=progress)
  }  
  class(out) <- c("mr.mash", "list")
  
  cat("Done!\n")
  toc <- Sys.time()
  cat("mr.mash successfully executed in", difftime(toc,tic, units="mins"),
      "minutes!\n")
  
  return(out)
}
