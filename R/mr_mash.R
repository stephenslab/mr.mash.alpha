#' @title  Multiple Regression with Multivariate Adaptive Shrinkage.
#' @details Performs multivariate multiple regression with mixture-of-normals prior.
#' 
#' @param Y an NxR matrix of responses.
#' @param X an NxP matrix of covariates.
#' @param V an RxR residual covariance matrix.
#' @param S0 a list of length K containing the desired RxR prior covariance matrices 
#' on the regression coefficients.
#' @param w0 a K-vector with prior mixture weights, each associated with the 
#' respective covariance matrix in \code{S0}.   
#' @param mu_init a PxR matrix of initial estimates for the regression coefficients.
#' @param tol convergence tolerance.
#' @param max_iter maximum number of iterations for the optimization algorithm.
#' @param update_w0 if TRUE, prior weights are updated.
#' @param compute_ELBO if TRUE, ELBO is computed.
#' 
#' @return a mr.mash fit, which is a list with some or all of the following elements\cr
#' \item{mu1}{a PxR matrix of posterior means for the regression coeffcients}
#' \item{S1}{a RxRxP array of posterior covariances for the regression coeffcients}
#' \item{w1}{a PxK matrix of posterior assignment probabilities to the mixture components}
#' \item{intercept}{an R-vector with the estimated intercepts}
#' \item{ELBO}{the Evidence Lower Bound at convergence} 
#'
#' @export
mr.mash <- function(Y, X, V, S0, w0, mu_init = matrix(0, nrow=ncol(X), ncol=ncol(Y)), 
                    tol=1e-8, max_iter=1e5, update_w0=T, compute_ELBO=T) {
  ###Center Y and X
  Y <- scale(Y, center=T, scale=F)
  X <- scale(X, center=T, scale=F)
 
  if(compute_ELBO){ 
    ###Compute inverse of V (needed for the ELBO)
    Vinv <- solve(V)
  }
  
  ###Initilize mu1, S1, w1, error, ELBO and iterator
  p       <- ncol(X)
  n       <- nrow(X)
  R       <- ncol(Y)
  K       <- length(S0)
  mu1_t   <- mu_init 
  S1_t    <- array(0, c(R, R, p))
  w1_t    <- matrix(0, nrow=p, ncol=K)
  err     <- matrix(Inf, nrow=p, ncol=R)
  t       <- 0
  if(compute_ELBO){
    ELBO    <- -Inf
  }
  
  ###Repeat the following until convergence
  cat("iter beta_max.diff ELBO_diff ELBO\n")
  while(any(err>tol)){
    ##Initialize ELBO parameters
    var_part_ERSS <- 0
    neg_KL <- 0
    
    if(compute_ELBO){
      ##Set last value of ELBO as ELBO0
      ELBO0 <- ELBO
    }
    
    ##Update iterator
    t <- t+1
    
    ##Exit loop if maximum number of iterations is reached
    if(t>max_iter){
      warning("Max number of iterations reached. Try increasing max_iter.")
      break
    }
    
    ##Compute expected residuals
    rbar <- Y - X%*%mu1_t
    
    ##Save current estimates.
    mu1_tminus1 <- mu1_t
    
    updates <- inner_loop(X=X, rbar=rbar, mu=mu_t, V=V, Vinv=Vinv, w0=w0, S0=S0, compute_ELBO=T)
    mu_t    <- updates$mu1
    S1_t    <- updates$S1
    w1_t    <- updates$w1
    rbar    <- updates$rbar
    
    ##Update w0 if requested
    if(update_w0){
      w0 <- update_weights(w1_t)
    }
    
    ##Compute distance in mu1 between two successive iterations
    err <- abs(mu1_t - mu1_tminus1)
    
    if(compute_ELBO){
      ##Compute ELBO
      var_part_ERSS <- updates$var_part_ERSS
      neg_KL <- updates$neg_KL
      ELBO <- compute_ELBO_fun(rbar=rbar, V=V, Vinv=Vinv, var_part_ERSS=var_part_ERSS, neg_KL=neg_KL)

      ##Print out useful info
      cat(sprintf("%4d %0.2e %0.2e %0.20e\n", t, max(err), ELBO - ELBO0, ELBO))
    } else {
      ##Print out useful info
      cat(sprintf("%4d %0.2e\n", t, max(err)))
    }
  }
  
  ###Compute intercept
  Ybar <- attr(Y, 'scaled:center')
  Xbar <- matrix(rep(attr(X, 'scaled:center'), each=ncol(mu1_t)), ncol=ncol(mu1_t), byrow=T)
  intercept <- Ybar - colSums(Xbar * mu1_t)
  
  if(compute_ELBO){
    ###Return the posterior assignment probabilities (w1), the posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), the intercept (intercept), and the Evidence Lower Bound (ELBO).    
    return(list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept, ELBO=ELBO))
  } else {
    ###Return the posterior assignment probabilities (w1), the posterior mean of the coefficients (mu1), and the posterior
    ###covariance of the coefficients (S1), and the intercept (intercept).    
    return(list(mu1=mu1_t, S1=S1_t, w1=w1_t, intercept=intercept))
  }
}
