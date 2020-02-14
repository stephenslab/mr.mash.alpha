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
#' @param scale_X if TRUE, X is centered and scaled. Scaling X allows a faster algorithm,
#' but the prior has a different interpretation.
#' @param verbose if TRUE, some information is printed to screen at each iteration.
#' 
#' @return a mr.mash fit, which is a list with some or all of the following elements\cr
#' \item{mu1}{a PxR matrix of posterior means for the regression coeffcients}
#' \item{S1}{a RxRxP array of posterior covariances for the regression coeffcients}
#' \item{w1}{a PxK matrix of posterior assignment probabilities to the mixture components}
#' \item{intercept}{an R-vector with the estimated intercepts}
#' \item{ELBO}{the Evidence Lower Bound at convergence}
#' 
#' @examples 
#' ###Set options
#' options(stringsAsFactors = F)
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
#'                rep(0, (p-2)*2)), byrow=T, ncol=2)
#'
#' ##Simulate X
#' X <- matrix(rnorm(n*p), nrow=n, ncol=p)
#' X <- scale(X, center=T, scale=F)
#'
#' ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
#' Y <- mr.mash.alpha:::sim_mvr(X, B, V)
#'
#' ###Specify the mixture weights and covariance matrices for the mixture-of-normals prior.
#' grid <- seq(1, 5)
#' S0mix <- mr.mash.alpha:::compute_cov_canonical(ncol(Y), singletons=T, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid, zeromat=T)
#' w0    <- rep(1/(length(S0mix)), length(S0mix))
#'
#' ###Estimate residual covariance
#' V_est <- cov(Y)
#'
#' ###Fit mr.mash
#' fit <- mr.mash(Y, X, V_est, S0mix, w0, tol=1e-8, update_w0=T, compute_ELBO=T, scale_X=T)
#'
#' @export
mr.mash <- function(Y, X, V, S0, w0, mu_init = matrix(0, nrow=ncol(X), ncol=ncol(Y)), 
                    tol=1e-8, max_iter=1e5, update_w0=T, compute_ELBO=T, scale_X=T,
                    verbose=T) {
  ###Center Y and either center and/or scale X
  Y <- scale(Y, center=T, scale=F)
  if(scale_X){
    X <- scale(X, center=T, scale=T)
  } else {
    X <- scale(X, center=T, scale=F)
  }
 
  ###Initilize mu1, S1, w1, error, ELBO and iterator
  p       <- ncol(X)
  n       <- nrow(X)
  R       <- ncol(Y)
  K       <- length(S0)
  mu1_t   <- mu_init 
  err     <- matrix(Inf, nrow=p, ncol=R)
  t       <- 0
  if(compute_ELBO){
    ELBO    <- -Inf
  }
  
  ###Precompute quantities
  if(scale_X){
    comps <- precompute_quants_scaled_X(n, V, S0)
  } else {
    comps <- precompute_quants_transformed_X(X, V, S0) 
  }
  
  if(compute_ELBO){ 
    ###Compute inverse of V (needed for the ELBO)
    Vinv <- chol2inv(comps$V_chol)
  }
  
  ###First iteration
  if(verbose){
    cat("iter beta_max.diff ELBO_diff ELBO\n")
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
  if(scale_X){
    ups   <- mr_mash_update_scaled_X(Y=Y, X=X, mu1_t=mu1_t, w1_t=NULL, V=V, Vinv=Vinv, ldetV=comps$ldetV, w0=w0, S0=S0, 
                                     S=comps$S, S1=comps$S1, SplusS0_chol=comps$SplusS0_chol, S_chol=comps$S_chol,
                                     ldetSplusS0_chol=comps$ldetSplusS0_chol, ldetS_chol=comps$ldetS_chol, 
                                     update_w0=update_w0, compute_ELBO=compute_ELBO)
  } else {
    ups   <- mr_mash_update(Y=Y, X=X, mu1_t=mu1_t, w1_t=NULL, V=V, Vinv=Vinv, ldetV=comps$ldetV, w0=w0, S0=S0, 
                            xtx=comps$xtx, V_chol=comps$V_chol, U0=comps$U0, d=comps$d, Q=comps$Q,
                            update_w0=update_w0, compute_ELBO=compute_ELBO)  
  }
  mu1_t <- ups$mu1_t
  S1_t  <- ups$S1_t
  w1_t  <- ups$w1_t
  if(compute_ELBO){
    ELBO  <- ups$ELBO
  }
  
  if(verbose){
    if(compute_ELBO){
      ##Print out useful info
      cat(sprintf("%4d %0.2e %0.2e %0.20e\n", t, max(err), ELBO - ELBO0, ELBO))
    } else {
      ##Print out useful info
      cat(sprintf("%4d %0.2e\n", t, max(err)))
    }
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
    if(scale_X){
      ups   <- mr_mash_update_scaled_X(Y=Y, X=X, mu1_t=mu1_t, w1_t=w1_t, V=V, Vinv=Vinv, ldetV=comps$ldetV, w0=w0, S0=S0, 
                                       S=comps$S, S1=comps$S1, SplusS0_chol=comps$SplusS0_chol, S_chol=comps$S_chol,
                                       ldetSplusS0_chol=comps$ldetSplusS0_chol, ldetS_chol=comps$ldetS_chol,
                                       update_w0=update_w0, compute_ELBO=compute_ELBO)
    } else {
      ups   <- mr_mash_update(Y=Y, X=X, mu1_t=mu1_t, w1_t=w1_t, V=V, Vinv=Vinv, ldetV=comps$ldetV, w0=w0, S0=S0, 
                              xtx=comps$xtx, V_chol=comps$V_chol, U0=comps$U0, d=comps$d, Q=comps$Q,
                              update_w0=update_w0, compute_ELBO=compute_ELBO)
    }
    mu1_t <- ups$mu1_t
    S1_t  <- ups$S1_t
    w1_t  <- ups$w1_t
    if(compute_ELBO){
      ELBO  <- ups$ELBO
    }
    
    ##Compute distance in mu1 between two successive iterations
    err <- abs(mu1_t - mu1_tminus1)
    
    if(verbose){
      if(compute_ELBO){
        ##Print out useful info
        cat(sprintf("%4d %0.2e %0.2e %0.20e\n", t, max(err), ELBO - ELBO0, ELBO))
      } else {
        ##Print out useful info
        cat(sprintf("%4d %0.2e\n", t, max(err)))
      }
    }
  }
  
  if(scale_X){
    ###Rescale posterior means and covariance of coefficients
    SX <- matrix(rep(attr(X, 'scaled:scale'), each=ncol(mu1_t)), ncol=ncol(mu1_t), byrow=T)
    mu1_t <- mu1_t/SX
    for(j in 1:dim(S1_t)[3]){
      S1_t[, , j] <- S1_t[, , j]/((attr(X, 'scaled:scale')[j])^2)
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
