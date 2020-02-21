#' @title  Multiple Regression with Multivariate Adaptive Shrinkage.
#' @details Performs multivariate multiple regression with mixture-of-normals prior.
#' 
#' @param Y an n x r matrix of responses.
#' @param X an n x p matrix of covariates.
#' @param V an r x r residual covariance matrix.
#' @param S0 a list of length K containing the desired r x r prior covariance matrices 
#' on the regression coefficients.
#' @param w0 a K-vector with prior mixture weights, each associated with the 
#' respective covariance matrix in \code{S0}.   
#' @param mu_init a p x r matrix of initial estimates for the regression coefficients.
#' @param tol convergence tolerance.
#' @param max_iter maximum number of iterations for the optimization algorithm.
#' @param update_w0 if TRUE, prior weights are updated.
#' @param compute_ELBO if TRUE, ELBO is computed.
#' @param standardize if TRUE, X is standardized using the sample means and sample standard deviations. 
#' Standardizing X allows a faster implementation, but the prior has a different interpretation.
#' Coefficients and covariances are returned on the original scale.
#' @param verbose if TRUE, some information is printed to screen at each iteration.
#' 
#' @return a mr.mash fit, which is a list with some or all of the following elements\cr
#' \item{mu1}{a p x r matrix of posterior means for the regression coeffcients.}
#' \item{S1}{a r x r x p array of posterior covariances for the regression coeffcients.}
#' \item{w1}{a p x K matrix of posterior assignment probabilities to the mixture components.}
#' \item{intercept}{an r-vector with the estimated intercepts.}
#' \item{fitted}{an n x r matrix of fitted values.}
#' \item{ELBO}{the Evidence Lower Bound at convergence.}
#' \item{progress}{A data.frame including information regarding convergence criteria at each iteration.}
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
#' ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
#' Y <- mr.mash.alpha:::sim_mvr(X, B, V)
#'
#' ###Specify the mixture weights and covariance matrices for the mixture-of-normals prior.
#' grid <- seq(1, 5)
#' S0mix <- mr.mash.alpha:::compute_cov_canonical(ncol(Y), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid, zeromat=TRUE)
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
#' fit <- mr.mash(Ytrain, Xtrain, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE, compute_ELBO=TRUE, standardize=TRUE)
#'
#' # Compare the "fitted" values of Y against the true Y in the training set.
#' plot(fit$fitted,Ytrain,pch = 20,col = "darkblue",xlab = "true",ylab = "fitted")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted")
#'
#' # Predict the multivariate outcomes in the test set using the fitted model.
#' Ytest_est <- predict(fit,Xtest)
#' plot(Ytest_est,Ytest,pch = 20,col = "darkblue",xlab = "true",ylab = "predicted")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted")
#' 
#' @export
mr.mash <- function(Y, X, V, S0, w0, mu_init=NULL, 
                    tol=1e-8, max_iter=1e5, update_w0=TRUE, compute_ELBO=TRUE, standardize=TRUE,
                    verbose=TRUE) {

  tic <- Sys.time()
  cat("Processing the inputs... ")
  
  ###Check that the inputs are in the correct format
  if(!is.matrix(Y)){
    stop("Y must be a matrix.")
  }
  if(!is.matrix(X)){
    stop("X must be a matrix.")
  }
  if (any(is.na(Y))) {
    stop("Y must not contain missing values.")
  }
  if (any(is.na(X))) {
    stop("X must not contain missing values.")
  }
  if(!is.matrix(V) && (isSymmetric(V))){
    stop("V must be a symmetric matrix.")
  }
  if(!is.list(S0)){
    stop("S0 must be a list.")
  }
  if(!is.vector(w0)){
    stop("w0 must be a vector.")
  }
  if(length(S0)!=length(w0)){
    stop("S0 and w0 must have the same length")
  }

  ###Center Y and either center and/or scale X
  Y <- scale(Y, center=TRUE, scale=FALSE)
  if(standardize){
    X <- scale(X, center=TRUE, scale=TRUE)
  } else {
    X <- scale(X, center=TRUE, scale=FALSE)
  }
 
  ###Initilize mu1, S1, w1, error, ELBO, iterator, and progress
  if(is.null(mu_init)){
    mu_init <- matrix(0, nrow=ncol(X), ncol=ncol(Y))
  }
  p        <- ncol(X)
  n        <- nrow(X)
  R        <- ncol(Y)
  K        <- length(S0)
  mu1_t    <- mu_init 
  err      <- matrix(Inf, nrow=p, ncol=R)
  t        <- 0
  progress <- data.frame() 
  if(compute_ELBO){
    ELBO    <- -Inf
  }
  
  ###Precompute quantities
  if(standardize){
    comps <- precompute_quants_scaled_X(n, V, S0)
  } else {
    comps <- precompute_quants_centered_X(X, V, S0) 
  }
  
  if(compute_ELBO){ 
    ###Compute inverse of V (needed for the ELBO)
    Vinv <- chol2inv(comps$V_chol)
    ldetV <- chol2ldet(comps$V_chol)
  } else {
    Vinv <- NULL
    ldetV <- NULL
  }
  
  cat("Done!\n")
  
  ###First iteration
  if(verbose){
    cat("Fitting the VEM algorithm... \n")
    if(compute_ELBO){
      cat(" iter    mu1_max.diff     ELBO_diff               ELBO\n")
    } else {
      cat(" iter    mu1_max.diff\n")
    }
  } else {
    cat("Fitting the VEM algorithm... ")
  }
  ##Save current estimates.
  mu1_tminus1 <- mu1_t   
  
  ##Update iterator
  t <- t+1
  
  if(compute_ELBO){
    ##Set last value of ELBO as ELBO0
    ELBO0 <- ELBO
  }
  
  ###Update variational parameters
  ups   <- mr_mash_update_general(Y=Y, X=X, mu1_t=mu1_t, w1_t=NULL, V=V, Vinv=Vinv, ldetV=ldetV, w0=w0, S0=S0, 
                                  precomp_quants=comps, update_w0=update_w0, compute_ELBO=compute_ELBO,
                                  standardize=standardize)
  mu1_t <- ups$mu1_t
  S1_t  <- ups$S1_t
  w1_t  <- ups$w1_t
  if(compute_ELBO){
    ELBO  <- ups$ELBO
  }
  
  if(compute_ELBO){
    if(verbose){
      ##Print out useful info
      cat(sprintf("%4d      %9.2e      %9.2e      %0.20e\n", t, max(err), ELBO - ELBO0, ELBO))
    }
    ##Update progress data.frame 
    progress <- rbind(progress, c(t, max(err), ELBO - ELBO0, ELBO))
  } else {
    if(verbose){
      ##Print out useful info
      cat(sprintf("%4d      %9.2e\n", t, max(err)))
    }
    ##Update progress data.frame 
    progress <- rbind(progress, c(t, max(err)))
  }
  
  ###Repeat the following until convergence
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
    
    if(compute_ELBO){
      ##Set last value of ELBO as ELBO0
      ELBO0 <- ELBO
    }
    
    ###Update model parameters and variational parameters
    ups   <- mr_mash_update_general(Y=Y, X=X, mu1_t=mu1_t, w1_t=w1_t, V=V, Vinv=Vinv, ldetV=ldetV, w0=w0, S0=S0, 
                                    precomp_quants=comps, update_w0=update_w0, compute_ELBO=compute_ELBO,
                                    standardize=standardize)
    mu1_t <- ups$mu1_t
    S1_t  <- ups$S1_t
    w1_t  <- ups$w1_t
    if(compute_ELBO){
      ELBO  <- ups$ELBO
    }
    
    ##Compute distance in mu1 between two successive iterations
    err <- abs(mu1_t - mu1_tminus1)
    
    if(compute_ELBO){
      if(verbose){
        ##Print out useful info
        cat(sprintf("%4d      %9.2e      %9.2e      %0.20e\n", t, max(err), ELBO - ELBO0, ELBO))
      }
      ##Update progress data.frame 
      progress <- rbind(progress, c(t, max(err), ELBO - ELBO0, ELBO))
    } else {
      if(verbose){
        ##Print out useful info
        cat(sprintf("%4d      %9.2e\n", t, max(err)))
      }
      ##Update progress data.frame 
      progress <- rbind(progress, c(t, max(err)))
    }
  }
  
  cat("Done!\n")
  cat("Processing the output... ")
  
  ###Compute fitted values
  fitted_vals <- X%*%mu1_t + matrix(rep(attr(Y,"scaled:center"), each=nrow(Y)), ncol=ncol(Y))
  attr(fitted_vals, "scaled:center") <- NULL
  attr(fitted_vals, "scaled:scale") <- NULL

  if(standardize){
    ###Rescale posterior means and covariance of coefficients
    SX <- matrix(rep(attr(X, 'scaled:scale'), each=ncol(mu1_t)), ncol=ncol(mu1_t), byrow=TRUE)
    mu1_t <- mu1_t/SX
    for(j in 1:dim(S1_t)[3]){
      S1_t[, , j] <- S1_t[, , j]/((attr(X, 'scaled:scale')[j])^2)
    }
  }

  ###Compute intercept
  intercept <- attr(Y,"scaled:center") - attr(X,"scaled:center") %*% mu1_t
  intercept <- drop(intercept)
  
  if(compute_ELBO){
    colnames(progress) <- c("iter", "mu1_max.diff", "ELBO_diff", "ELBO")
    
    ###Return the posterior assignment probabilities (w1), the posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the intercept (intercept), the fitted values (fitted), the Evidence Lower Bound (ELBO), 
    ###and the progress data frame (progress).
    out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept, fitted=fitted_vals, ELBO=ELBO, progress=progress)
  } else {
    colnames(progress) <- c("iter", "mu1_max.diff")
    
    ###Return the posterior assignment probabilities (w1), the posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the intercept (intercept), the fitted values (fitted), and the progress data frame (progress).    
    out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept, fitted=fitted_vals, progress=progress)
  }
  
  class(out) <- c("mr.mash", "list")
  
  cat("Done!\n")
  toc <- Sys.time()
  cat("mr.mash successfully executed in", difftime(toc,tic, units="mins"), "minutes!\n")
  
  return(out)
}
