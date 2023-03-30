#' @title  Multiple Regression with Multivariate Adaptive Shrinkage
#'  from summary data.
#' 
#' @description Performs multivariate multiple regression with
#'   mixture-of-normals prior.
#' 
#' @param Bhat p x r matrix of regression coefficients from univariate
#'  simple linear regression.
#' 
#' @param Shat p x r matrix of standard errors of the regression coefficients 
#'  from univariate simple linear regression.
#'  
#' @param Z p x r matrix of Z-scores from univariate
#'  simple linear regression.
#'  
#' @param R p x p correlation matrix among the variables.
#' 
#' @param covY r x r covariance matrix across responses.
#' 
#' @param n scalar indicating the sample size.
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
#'   mean regression coefficients. These should be on the same scale as
#'   the X provided. If \code{standardize=TRUE}, mu1_init will be scaled
#'   appropriately after standardizing X.
#'   
#' @param convergence_criterion Criterion to use for convergence check.
#' 
#' @param tol Convergence tolerance.
#' 
#' @param max_iter Maximum number of iterations for the optimization
#'   algorithm.
#' 
#' @param update_w0 If \code{TRUE}, prior weights are updated.
#' 
#' @param update_w0_method Method to update prior weights. Only EM is
#'   currently supported.
#' 
#' @param w0_threshold Drop mixture components with weight less than this value.
#'   Components are dropped at each iteration after 15 initial iterations.
#'   This is done to prevent from dropping some poetentially important 
#'   components prematurely.
#' 
#' @param update_V if \code{TRUE}, residual covariance is updated.
#' 
#' @param update_V_method Method to update residual covariance. So far,
#'   "full" and "diagonal" are supported. If \code{update_V=TRUE} and V 
#'   is not provided by the user, this option will determine how V is 
#'   computed (and fixed) internally from \code{mu1_init}.
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
#' @param ca_update_order The order with which coordinates are
#'   updated.  So far, "consecutive", "decreasing_logBF",
#'   "increasing_logBF", "random" are supported.
#'   
#' @param X_colmeans a p-vector of variable means.
#' 
#' @param Y_colmeans a r-vector of response means.
#' 
#' @param check_R If \code{TRUE}, R is checked to be positive semidefinite.
#' 
#' @param R_tol tolerance to declare positive semi-definiteness of R.
#'   
#' @param nthreads Number of RcppParallel threads to use for the
#'   updates. When \code{nthreads} is \code{NA}, the default number of
#'   threads is used; see
#'   \code{\link[RcppParallel]{defaultNumThreads}}. This setting is
#'   ignored when \code{version = "R"}.
#' 
#' @return A mr.mash.rss fit, stored as a list with some or all of the
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
#'   intercept, if \code{X_colmeans} and \code{Y_colmeans} are provided.
#'   Otherwise, \code{NA} is output.}
#' 
#' \item{ELBO}{Evidence Lower Bound (ELBO) at last iteration.}
#' 
#' \item{progress}{A data frame including information regarding
#'   convergence criteria at each iteration.}
#'   
#' \item{converged}{\code{TRUE} or \code{FALSE}, indicating whether
#'   the optimization algorithm converged to a solution within the chosen tolerance
#'   level.}
#'   
#' \item{Y}{n x r matrix of responses at last iteration (only relevant when missing values
#'   are present in the input Y).}
#'  
#' @examples 
#' ###Set seed
#' set.seed(123)
#'
#' ###Simulate X and Y
#' ##Set parameters
#' n  <- 1000
#' p <- 100
#' p_causal <- 20
#' r <- 5
#'
#' ###Simulate data
#' out <- simulate_mr_mash_data(n, p, p_causal, r, pve=0.5, B_cor=1,
#'                              B_scale=1, X_cor=0, X_scale=1, V_cor=0)
#' 
#' ###Split the data in training and test sets
#' Ytrain <- out$Y[-c(1:200), ]
#' Xtrain <- out$X[-c(1:200), ]
#' Ytest <- out$Y[c(1:200), ]
#' Xtest <- out$X[c(1:200), ]
#' 
#' ###Specify the covariance matrices for the mixture-of-normals prior.
#' univ_sumstats <- compute_univariate_sumstats(Xtrain, Ytrain,
#'                    standardize=TRUE, standardize.response=FALSE)
#' grid <- autoselect.mixsd(univ_sumstats, mult=sqrt(2))^2
#' S0 <- compute_canonical_covs(ncol(Ytrain), singletons=TRUE,
#'                              hetgrid=c(0, 0.25, 0.5, 0.75, 1))
#' S0 <- expand_covs(S0, grid, zeromat=TRUE)
#'
#' ###Fit mr.mash
#' covY <- cov(Ytrain)
#' corX <- cor(Xtrain)
#' n_train <- nrow(Ytrain)
#' fit <- mr.mash.rss(Bhat=univ_sumstats$Bhat, Shat=univ_sumstats$Shat, S0=S0, 
#'                    covY=covY, R=corX, n=n_train, V=covY, update_V=TRUE,
#'                    X_colmeans=colMeans(Xtrain), Y_colmeans=colMeans(Ytrain))
#'
#' # Predict the multivariate outcomes in the test set using the fitted model.
#' Ytest_est <- predict(fit,Xtest)
#' plot(Ytest_est,Ytest,pch = 20,col = "darkblue",xlab = "true",
#'      ylab = "predicted")
#' abline(a = 0,b = 1,col = "magenta",lty = "dotted")
#'
#' @importFrom stats cov
#' @importFrom Rfast is.symmetric
#' @importFrom RcppParallel defaultNumThreads
#' @importFrom RcppParallel setThreadOptions
#'
#' @export
#' 
mr.mash.rss <- function(Bhat, Shat, Z, R, covY, n, S0, w0=rep(1/(length(S0)), length(S0)), V=NULL, 
                        mu1_init=NULL, tol=1e-4, convergence_criterion=c("mu1", "ELBO"),
                        max_iter=5000, update_w0=TRUE, update_w0_method="EM", 
                        w0_threshold=0, compute_ELBO=TRUE, standardize=TRUE, verbose=TRUE,
                        update_V=FALSE, update_V_method=c("full", "diagonal"), version=c("Rcpp", "R"), e=1e-8,
                        ca_update_order=c("consecutive", "decreasing_logBF", "increasing_logBF", "random"),
                        X_colmeans=NULL, Y_colmeans=NULL, check_R=TRUE, R_tol=1e-08,
                        nthreads=as.integer(NA)) {
  
  if(verbose){
    tic <- Sys.time()
    cat("Processing the inputs... ")
  }

  # CHECK AND PROCESS INPUTS
  # ------------------------
  ###Select method to check for convergence (if not specified by user, mu1
  ###will be used)
  convergence_criterion <- match.arg(convergence_criterion)
  
  ###Select method to update the weights (if not specified by user, EM
  ###will be used)
  update_w0_method <- match.arg(update_w0_method)
  
  ###Select method to update the residual covariance (if not specified by user, full
  ###will be used)
  update_V_method <- match.arg(update_V_method)
  
  ###Select version of the inner loop (if not specified by user, Rcpp
  ###will be used)
  version <- match.arg(version)
  
  ###Select ordering of the coordinate ascent updates (if not specified by user,
  ###consecutive will be used
  ca_update_order <- match.arg(ca_update_order)
  
  ###Initialize the RcppParallel multithreading using a pre-specified number
  ###of threads, or using the default number of threads when nthreads is NA.
  if(version=="Rcpp"){
    if (is.na(nthreads)) {
      setThreadOptions()
      nthreads <- defaultNumThreads()
    } else
      setThreadOptions(numThreads = nthreads)
  }
  
  ###Check that the inputs are in the correct format
  if (sum(c(missing(Z), missing(Bhat) || missing(Shat))) != 1)
    stop("Please provide either Z or (Bhat, Shat), but not both")
  
  if(missing(Z)){
    if(!is.matrix(Bhat))
      stop("Bhat must be a matrix.")
    if(!is.matrix(Shat))
      stop("Shat must be a matrix.")
    if(any(is.na(Bhat)) || any(is.na(Shat)))
      stop("Bhat, Shat must not contain missing values.")
    if(any(Shat <= 0))
      stop("Shat cannot have zero or negative elements.")
  } else {
    if(!is.matrix(Z))
      stop("Z must be a matrix.")
    if(any(is.na(Z)))
      stop("Z must not contain missing values.")
  }
  
  if(!is.null(V)){
    if(!is.matrix(V) || !is.symmetric(V))
      stop("V must be a symmetric matrix.")
  }
  if(!missing(covY)){
    if(!is.matrix(covY) || !is.symmetric(covY))
      stop("covY must be a symmetric matrix.")
  }
  if(!is.matrix(R) || !is.symmetric(R))
    stop("R must be a symmetric matrix.")
  if(!is.list(S0))
    stop("S0 must be a list.")
  if(!is.vector(w0))
    stop("w0 must be a vector.")
  if(length(w0)<2)
    stop("At least 2 mixture components must be present.")
  if(abs(sum(w0) - 1) > 1e-10)
    stop("Elements of w0 must sum to 1.")
  if(length(S0)!=length(w0))
    stop("S0 and w0 must have the same length.")
  if(!missing(mu1_init) && !is.matrix(mu1_init))
    stop("mu1_init must be a matrix.")
  if(convergence_criterion=="ELBO" && !compute_ELBO)
    stop("ELBO needs to be computed with convergence_criterion=\"ELBO\".")

  # PRE-PROCESSING STEPS
  # --------------------
  
  ###Compute Z scores
  if(missing(Z)){
    Z <- Bhat/Shat
  }
  
  Z[is.na(Z)] <- 0
  
  ###Compute pve-adjusted Z scores, if n is provided
  if(!missing(n)) {
    adj <- (n-1)/(Z^2 + n - 2)
    Z   <- sqrt(adj) * Z
  }
  
  ###Obtain dimensions and store dimensions names of the inputs
  p <- nrow(Z)
  r <- ncol(Z)
  K <- length(S0)
  Z_colnames <- colnames(Z)
  Z_rownames <- rownames(Z)
  
  ###If covariance of Y and standard errors are provided,
  ###the effects are on the *original scale*.
  if(!missing(Shat) & !missing(covY)){
    XtXdiag <- rowMeans(matrix(diag(covY), nrow=p, ncol=r, byrow=TRUE) * adj/(Shat^2))
    XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
    XtX <- (XtX + t(XtX))/2
    XtY <- Z * sqrt(adj) * matrix(diag(covY), nrow=p, ncol=r, byrow=TRUE) / Shat
  } else {
    ###The effects are on the *standardized* X, y scale.
    XtX <- R*(n-1)
    XtY <- Z*sqrt(n-1)
    covY <- cov2cor(V)
  }
  
  YtY <- covY*(n-1)
  
  ###Check whether XtX is positive semidefinite
  if(check_R){
    semi_pd <- check_semi_pd(XtX, R_tol)
    if (!semi_pd$status)
      stop("XtX is not a positive semidefinite matrix")
  }
  
  ###Adjust XtX and XtY if X is standardized 
  if(standardize){
    dXtX <- diag(XtX)
    sx <- sqrt(dXtX/(n-1))
    sx[sx == 0] <- 1
    XtX <- t((1/sx) * XtX) / sx
    XtY <-  XtY / sx
  }
  
  ###Add number to diagonal elements of the prior matrices (improves
  ###numerical stability)
  S0 <- lapply(S0, makePD, e=e)
  
  ###Initialize regression coefficients to 0 if not provided
  if(is.null(mu1_init)){
    mu1_init <- matrix(0, nrow=p, ncol=r)
  }
  
  ###Scale mu1_init, if X is standardized 
  if(standardize)
    mu1_init <- mu1_init*sx 
  
  ###Initilize mu1, S1, w1, delta_mu1, delta_ELBO, delta_conv, ELBO, iterator, progress
  mu1_t <- mu1_init 
  delta_mu1 <- matrix(Inf, nrow=p, ncol=r)
  delta_ELBO <- Inf
  if(convergence_criterion=="mu1")
    delta_conv <- max(delta_mu1)
  else if(convergence_criterion=="ELBO")
    delta_conv <- delta_ELBO
  ELBO <- -Inf
  t <- 0
  progress <- as.data.frame(matrix(as.numeric(NA), nrow=max_iter, ncol=3))
  colnames(progress) <- c("iter", "timing", "mu1_max.diff")
  if(compute_ELBO){
    progress$ELBO_diff <- as.numeric(NA)
    progress$ELBO <- as.numeric(NA)
  }

  ###Compute V, if not provided by the user
  if(is.null(V)){
    # How to do so with sumstats??
    
    if(update_V_method=="diagonal")
      V <- diag(diag(V))
  }

  ###Set eps
  eps <- .Machine$double.eps
  
  ###Precompute quantities
  comps <- precompute_quants(n, V, S0, standardize, version)
  if(!standardize){
    comps$xtx <- diag(XtX)
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
    update_order <- compute_rank_variables_BFmix_rss(n, XtY, V, Vinv, w0, S0, comps, standardize, version,
                                                     decreasing=TRUE, eps, nthreads)
  } else if(ca_update_order=="increasing_logBF"){
    update_order <- compute_rank_variables_BFmix_rss(n, XtY, V, Vinv, w0, S0, comps, standardize, version,
                                                     decreasing=FALSE, eps, nthreads)
  } else if(ca_update_order=="random")
    update_order <- sample(x=1:p, size=p)
  
  if(verbose)
    cat("Done!\n")

  # MAIN LOOP
  # ---------
  if(verbose){
    if(version=="Rcpp" && nthreads>1)
      cat(sprintf("Fitting the optimization algorithm using %d RcppParallel threads... \n", nthreads))
    else
      cat("Fitting the optimization algorithm... \n")
    cat(" iter    mu1_max.diff")
    if(compute_ELBO)
      cat("     ELBO_diff               ELBO\n")
    else 
      cat("\n")
  }
  
  ###Repeat the following until convergence, or until maximum number
  ###of iterations is reached.
  while(delta_conv>tol){
    
    ##Start timing
    time1 <- proc.time()
      
    ##Save current estimates.
    mu1_old <- mu1_t   
    
    ##Set last value of ELBO as ELBO_old
    ELBO_old <- ELBO
    
    ##Update iterator
    t <- t+1
    
    ##Exit loop if maximum number of iterations is reached
    if(t>max_iter){
      warning("Max number of iterations reached. Try increasing max_iter.")
      break
    }

    # M-STEP
    # ------
    if(t > 1){
      ##Update V if requested
      if(update_V){
        V <- update_V_rss_fun(n, RbartRbar, var_part_ERSS)
        if(update_V_method=="diagonal")
          V <- diag(diag(V))

        #Recompute precomputed quantities after updating V
        comps <- precompute_quants(n, V, S0, standardize, version)
        if(!standardize)
          comps$xtx <- diag(XtX)
        if(compute_ELBO || !standardize)
          Vinv <- chol2inv(comps$V_chol)
        if(compute_ELBO)
          ldetV <- chol2ldet(comps$V_chol)
      }
      
      ##Update w0 if requested
      if(update_w0){
        w0 <- update_weights_em(w1_t)
        
        #Drop components with mixture weight <= w0_threshold
        if(t>15 && any(w0 < w0_threshold)){
          to_keep <- which(w0 >= w0_threshold)
          w0 <- w0[to_keep]
          w0 <- w0/sum(w0)
          S0 <- S0[to_keep]
          if(length(to_keep) > 1){
            comps <- filter_precomputed_quants(comps, to_keep, standardize, version)
          } else if(length(to_keep) == 1 & all((S0[[to_keep]] - (diag(nrow(S0[[to_keep]]))*e)) < eps)){ #null component is the only one left
            mu1_t <- matrix(0, nrow=p, ncol=r)
            S1_t <- array(0, c(r, r, p))
            w1_t <- matrix(1, nrow=p, ncol=1)
            warning("Only the null component is left. Estimated coefficients are set to 0.")
            break
          } else { #some other component is the only one left
            stop("Only one component (different from the null) left. Consider lowering w0_threshold.")
          }
	}
      }
    }
    
    # E-STEP
    # ------
    ###Update variational parameters
    ups <- mr_mash_update_general_rss(n=n, XtX=XtX, XtY=XtY, YtY=YtY, 
                                      mu1_t=mu1_t, V=V,
                                      Vinv=Vinv, ldetV=ldetV, w0=w0, S0=S0, 
                                      precomp_quants=comps,
                                      compute_ELBO=compute_ELBO,
                                      standardize=standardize,
                                      update_V=update_V, version=version, 
                                      update_order=update_order, eps=eps,
                                      nthreads=nthreads)
    mu1_t <- ups$mu1_t
    S1_t  <- ups$S1_t
    w1_t  <- ups$w1_t
    if(compute_ELBO)
      ELBO <- ups$ELBO
    if(update_V)
      var_part_ERSS <- ups$var_part_ERSS
    if(compute_ELBO || update_V)
      RbartRbar <- ups$RbartRbar

    ##End timing
    time2 <- proc.time()
    
    ##Compute difference in mu1 and ELBO between two successive iterations,
    ##and assign the requested criterion to delta_conv
    delta_mu1 <- abs(mu1_t - mu1_old)
    delta_ELBO <- ELBO - ELBO_old
    if(convergence_criterion=="mu1")
      delta_conv <- max(delta_mu1)
    else if(convergence_criterion=="ELBO")
      delta_conv <- delta_ELBO
    
    ##Update progress data.frame 
    progress[t, c(1:3)] <- c(t, time2["elapsed"] - time1["elapsed"], max(delta_mu1))
    if(compute_ELBO)
      progress[t, c(4, 5)] <- c(delta_ELBO, ELBO)
    
    if(verbose){
      ##Print out useful info
      cat(sprintf("%4d      %9.2e", t, max(delta_mu1)))
      if(compute_ELBO)
        cat(sprintf("      %9.2e      %0.20e\n", delta_ELBO, ELBO))
      else
        cat("\n")
    }
  }
  
  ###Record convergence status
  if(t>max_iter)
    converged <- FALSE
  else
    converged <- TRUE
  
  if(verbose){
    cat("Done!\n")
    cat("Processing the outputs... ")
  }

  # POST-PROCESSING STEPS
  # --------------------
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
  if(!is.null(X_colmeans) & !is.null(Y_colmeans)){
    intercept <- drop(Y_colmeans - X_colmeans %*% mu1_t)
    names(intercept) <- Z_colnames
  }
  
  ###Assign names to outputs dimensions
  S0_names <- names(S0)
  rownames(mu1_t) <- Z_rownames
  colnames(mu1_t) <- Z_colnames
  dimnames(S1_t)[[1]] <- Z_colnames
  dimnames(S1_t)[[2]] <- Z_colnames
  dimnames(S1_t)[[3]] <- Z_rownames
  S0 <- lapply(S0, function(x){rownames(x) <- colnames(x) <- Z_colnames; return(x)})
  rownames(w1_t) <- Z_rownames
  colnames(w1_t) <- S0_names
  names(w0) <- S0_names
  rownames(V) <- Z_colnames
  colnames(V) <- Z_colnames

  ###Remove unused rows of progress
  progress <- progress[rowSums(is.na(progress)) != ncol(progress), ]
  
  ###Return the posterior assignment probabilities (w1), the
  ###posterior mean of the coefficients (mu1), and the posterior
  ###covariance of the coefficients (S1), the residual covariance (V),
  ###the prior weights (w0), the intercept (intercept), the fitted values (fitted), 
  ###and the progress data frame (progress), the prior covariance (S0), convergence
  ###status, the covariance of the fitted values (G), the proportion of variance explained (pve),
  ###the Evidence Lower Bound (ELBO; if computed) and imputed responses (Y; if 
  ###missing values were present).
  out <- list(mu1=mu1_t, S1=S1_t, w1=w1_t, V=V, w0=w0, S0=simplify2array_custom(S0), 
              intercept=NA, progress=progress, converged=converged)
  if(compute_ELBO)
    ###Append ELBO to the output
    out$ELBO <- ELBO
  if(!is.null(X_colmeans) & !is.null(Y_colmeans))
    out$intercept <- intercept

  class(out) <- c("mr.mash.rss", "list")
  
  if(verbose){
    cat("Done!\n")
    toc <- Sys.time()
    cat("mr.mash.rss successfully executed in", difftime(toc, tic, units="mins"),
        "minutes!\n")
  }
  
  return(out)
}
