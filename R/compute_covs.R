#' @title Compute canonical covariance matrices.
#' @description Function to compute canonical covariance matrices scaled.
#' 
#' @param r number of responses.
#' 
#' @param singletons if \code{TRUE}, the response-specific effect matrices will be 
#'        included.
#'        
#' @param hetgrid scalar or numeric vector of positive correlation [0, 1] of the effects
#'        across responses.
#'
#' @return A list containing the canonical covariance matrices.
#'   
#' @export
compute_canonical_covs <- function(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1)){
  mats <- list()
  s_idx <- 0
  nms <- vector("character")
  ###Singleton matrices
  if(singletons) {
    for(i in 1:r) {
      mats[[i]] <- matrix(0, nrow=r, ncol=r)
      mats[[i]][i, i] <- 1
      nms[i] = paste0('singleton', i)
    }
    s_idx <- r
  }
  ###Heterogeneity matrices
  if(!is.null(hetgrid)) {
    for(j in 1:length(hetgrid)) {
      mats[[s_idx+j]] <- matrix(1, nrow=r, ncol=r)
      mats[[s_idx+j]][lower.tri(mats[[s_idx+j]], diag = FALSE)] <- hetgrid[j]
      mats[[s_idx+j]][upper.tri(mats[[s_idx+j]], diag = FALSE)] <- hetgrid[j]
      if(any(mats[[s_idx+j]] != diag(r)))
        nms[s_idx+j] <- paste0('shared', hetgrid[j])
      else
        nms[s_idx+j] <- 'independent'
    }
  }
  names(mats) <- nms
  return(mats)
}

#' @title Compute data-driven covariance matrices.
#' @description Function to compute data-driven covariance matrices from summary statistics
#'   using PCA, FLASH and the sample covariance. These matrices are de-noised using Extreme Deconvolution.
#' 
#' @param sumstats a list with two elements. 1 - Bhat, a numeric vector of regression coefficients.
#'   2 - Shat, a numeric vector of of standard erros for the regression coefficients.
#' 
#' @param subset_thresh scalar indicating the threshold for selecting the effects to be used for computing 
#'   the covariance matrices based on false local sign rate (lfsr) for a response-by-response ash analysis.  
#' 
#' @param n_pcs indicating the number of principal components to be selected.
#' 
#' @param flash_factors factors "default" to use \code{flashr} default function to initialize factors, currently \code{udv_si}. 
#' "nonneg" to implement a non-negative constraint on the factors
#' 
#' @param flash_remove_singleton whether or not factors corresponding to singleton matrices should be removed from output.
#' 
#' @param Gamma an r x r correlation matrix for the residuals; must be positive
#'   definite.
#'
#' @return A list containing the (de-noised) data-driven covariance matrices.
#' 
#' @importFrom mashr mash_set_data mash_1by1 get_significant_results cov_pca cov_flash cov_ed
#'   
#' @export
compute_data_driven_covs <- function(sumstats, subset_thresh=NULL, n_pcs=3, flash_factors=c("default", "nonneg"),
                                     flash_remove_singleton=FALSE,
                                     Gamma=diag(ncol(sumstats$Bhat))){
  
  flash_factors <- match.arg(flash_factors)
  
  ###Obtain strong effects
  data <- mash_set_data(sumstats$Bhat, sumstats$Shat, V=Gamma)
  if(!is.null(subset_thresh)){
    m_1by1 <- mash_1by1(data)
    subs <- get_significant_results(m_1by1, subset_thresh)
  } else {
    subs <- NULL
  }
  
  ##Compute data-driven matrices
  U_pca <- cov_pca(data=data, npc=n_pcs, subset=subs)
  U_flash <- cov_flash(data=data, factors=flash_factors, subset=subs, tag=NULL,
                       remove_singleton=flash_remove_singleton, output_model=NULL)
  U_emp <- cov_empirical(data=data, subset=subs)

  ##De-noise data-driven matrices via extreme deconvolution
  U_datadriven <- c(U_pca, U_flash, list(BB=U_emp))
  U_ed <- cov_ed(data, U_datadriven, subset=subs)
  
  return(U_ed)
}

#' @title Compute summary statistics from univariate simple linear regression.
#' @description Function to compute regression coefficients and their standard errors
#'   from univariate simple linear regression.
#' 
#' @param X n x p matrix of covariates.
#' 
#' @param Y n x r matrix of responses.
#'        
#' @param standardize If \code{TRUE}, X is "standardized" using the
#'   sample means and sample standard deviations.
#'        
#' @param standardize.response If \code{TRUE}, Y is "standardized" using the
#'   sample means and sample standard deviations.
#'   
#' @param mc.cores Number of cores to use. Parallelization is done over responses. 
#'   
#' @return A list with following elements:
#' 
#' \item{Bhat}{p x r matrix of the regression coeffcients.}
#'
#' \item{Shat}{p x r matrix of the standard errors for regression coeffcients.}
#' 
#' @importFrom parallel makeCluster parLapply stopCluster
#' 
#' @export
compute_univariate_sumstats <- function(X, Y, standardize=FALSE, standardize.response=FALSE, mc.cores=1){
  r <- ncol(Y)
  
  X <- scale_fast2(X, scale=standardize)$M 
  Y <- scale_fast2(Y, scale=standardize.response)$M
  
  linreg <- function(i, X, Y){
    p <- ncol(X)
    bhat <- rep(as.numeric(NA), p)
    shat <- rep(as.numeric(NA), p)
    
    for(j in 1:p){
      fit <- lm(Y[, i] ~ X[, j]-1)
      bhat[j] <- coef(fit)
      shat[j] <- summary(fit)$coefficients[1, 2]
    }
    
    return(list(bhat=bhat, shat=shat))
  }
  
  if(mc.cores>1){
    ###mclapply is a little faster but uses more memory
    # out <- mclapply(1:r, linreg, X, Y, mc.cores=mc.cores)
    
    cl <- makeCluster(mc.cores)
    out <- parLapply(cl, 1:r, linreg, X, Y)
    stopCluster(cl)
  } else {
    out <- lapply(1:r, linreg, X, Y)
  }
  
  return(list(Bhat=sapply(out,"[[","bhat"), Shat=sapply(out,"[[","shat")))
}

#' @title Compute a grid of standard deviations to scale the canonical covariance matrices.
#' @description Function to compute a grid of standard deviations from univariate 
#'   simple linear regression summary statistics
#' 
#' @param data a list with two elements. 1 - Bhat, a numeric vector of regression coefficients.
#'   2 - Shat, a numeric vector of of standard erros for the regression coefficients.
#' 
#' @param mult a scalar affecting how dense the resulting grid of standard deviations will be. 
#'   
#' @return A numeric vector of standard deviations.
#' 
#' @export
autoselect.mixsd <- function(data, mult=2){
  include <- !(data$Shat==0 | !is.finite(data$Shat) | is.na(data$Bhat))
  gmax <- grid_max(data$Bhat[include], data$Shat[include])
  gmin <- grid_min(data$Bhat[include], data$Shat[include])
  if (mult == 0) {
    return(c(0, gmax/2))
  }
  else {
    npoint = ceiling(log2(gmax/gmin)/log2(mult))
    return(mult^((-npoint):0) * gmax)
  }
}

#' @title Expand covariance matrices by a grid of variances for use in \code{mr.mash}.
#' @description Function to scale the covariance matrices by a grid of variances.
#' 
#' @param mats a list of covariance matrices.
#'        
#' @param grid scalar or numeric vector of variances of the effects.
#' 
#' @param zeromat if \code{TRUE}, the no-effect matrix will be included.
#'
#' @return A list containing the scaled covariance matrices.
#'   
#' @export
expand_covs <- function(mats, grid, zeromat=TRUE){
  r <- ncol(mats[[1]])
  mats <- lapply(mats, normalize_cov)
  U <- list()
  nms <- vector("character")
  t <- 0
  for(i in 1:length(mats)){
    for(j in 1:length(grid)){
      t <- t+1
      U[[t]] <- mats[[i]]*grid[j]
      nms[t] <- paste0(names(mats)[i], "_grid", j) 
    }
  }
  
  if(zeromat){
    zero_mat <- list(matrix(0, r, r))
    U <- c(zero_mat, U)
    nms <- c("null", nms)
  }
  
  names(U) <- nms
  
  return(U)
}

###Normalize a covariance matrix so its maximum diagonal element is 1.
normalize_cov <- function(U){
  if(max(diag(U))!=0){
    U = U/max(diag(U))
  }
  return(U)
}

###Compute the minimum value for the grid
grid_min = function(Bhat,Shat){
  min(Shat)/10
}

###Compute the maximum value for the grid
grid_max = function(Bhat,Shat){
  if (all(Bhat^2 <= Shat^2)) {
    8 * grid_min(Bhat,Shat) # the unusual case where we don't need much grid
  }  else {
    2 * sqrt(max(Bhat^2 - Shat^2))
  }
}

###Compute empirical covariance
cov_empirical <- function(data, subset=NULL){
  if(is.null(subset)) 
    subset <- 1:nrow(data$Bhat)
  B_center <- apply(data$Bhat[subset, ], 2, function(x) x - mean(x))
  U_emp <- crossprod(B_center) / nrow(B_center)
  
  return(U_emp)
}





#################################################
##N.B. functions below deprecated and left only## 
##for backwards compatibility                  ##
#################################################

###Function to compute canonical covariance matrices
create_cov_canonical <- function(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1)){
  mats <- list()
  
  ###Singleton matrices
  if((singletons)){
    for(i in 1:r){
      mats[[i]] <- matrix(0, nrow=r, ncol=r)
      mats[[i]][i, i] <- 1
    }
    
    ###Heterogeneity matrices
    if(!is.null(hetgrid)){
      for(j in 1:length(hetgrid)){
        mats[[r+j]] <- matrix(1, nrow=r, ncol=r)
        mats[[r+j]][lower.tri(mats[[r+j]], diag = FALSE)] <- hetgrid[j]
        mats[[r+j]][upper.tri(mats[[r+j]], diag = FALSE)] <- hetgrid[j]
      }
    }
  } else {
    ###Heterogeneity matrices
    if(!is.null(hetgrid)){
      for(j in 1:length(hetgrid)){
        mats[[j]] <- matrix(1, nrow=r, ncol=r)
        mats[[j]][lower.tri(mats[[j]], diag = FALSE)] <- hetgrid[j]
        mats[[j]][upper.tri(mats[[j]], diag = FALSE)] <- hetgrid[j]
      }
    }
  }
  return(mats)
}

###Compute canonical covariance matrices scaled by a grid 
compute_cov_canonical <- function(ntraits, singletons, hetgrid, grid, zeromat=TRUE){
  S <- create_cov_canonical(ntraits, singletons, hetgrid)
  U <- list()
  t <- 0
  for(i in 1:length(S)){
    for(j in 1:length(grid)){
      t <- t+1
      U[[t]] <- S[[i]]*grid[j]
    }
  }
  
  names(U) <- paste0("S0_", seq(1, length(U)))
  
  if(zeromat){
    zero_mat <- matrix(0, ntraits, ntraits)
    #zero_mat[upper.tri(zero_mat)] <- 1e-10
    #zero_mat[lower.tri(zero_mat)] <- 1e-10
    U[[paste0("S0_", length(U)+1)]] <- zero_mat
  }
  
  return(U)
}
###Compute univariate summary statistics
get_univariate_sumstats <- function(X, Y, standardize=FALSE, standardize.response=FALSE){
  r <- ncol(Y)
  p <- ncol(X)
  B <- matrix(as.numeric(NA), nrow=p, ncol=r)
  S <- matrix(as.numeric(NA), nrow=p, ncol=r)
  
  X <- scale(X, center=TRUE, scale=standardize) 
  Y <- scale(Y, center=TRUE, scale=standardize.response)
  
  for(i in 1:r){
    for(j in 1:p){
      fit <- lm(Y[, i] ~ X[, j]-1)
      B[j, i] <- coef(fit)
      S[j, i] <- summary(fit)$coefficients[1, 2]
    }
  }
  
  return(list(Bhat=B, Shat=S))
}

