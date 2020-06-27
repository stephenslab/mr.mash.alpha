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
#'   the covariance matrices based on false local sign rate (lfsr).  
#' 
#' @param n_pcs indicating the number of principal components to be selected.
#' 
#' @param non_canonical ???
#'
#' @return A list containing the (de-noised) data-driven covariance matrices.
#' 
#' @importFrom mashr mash_set_data mash_1by1 get_significant_results cov_pca cov_ed
#'   
#' @export
compute_data_driven_covs <- function(sumstats, subset_thresh=NULL, n_pcs=3, non_canonical=TRUE){
  ###Obtain string effects
  data <- mash_set_data(sumstats$Bhat, sumstats$Shat)
  if(!is.null(subset_thresh)){
    m_1by1 <- mash_1by1(data)
    subs <- get_significant_results(m_1by1, subset_thresh)
  } else {
    subs <- NULL
  }
  
  ##Compute data-driven matrices
  U_pca <- cov_pca(data=data, npc=n_pcs, subset=subs)
  U_flash <- cov_flash(data=data, subset=subs, non_canonical=non_canonical, save_model=NULL)
  if(is.null(subs)){
    subs = 1:mashr:::n_effects(data)
  }
  B_center <- apply(data$Bhat[subs, ], 2, function(x) x - mean(x))
  U_BB <- t(B_center) %*% B_center / nrow(B_center)
  
  ##De-noise data-driven matrices via extreme deconvolution
  U_datadriven <- c(U_pca, U_flash, list(BB=U_BB))
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
#' @return A list with following elements:
#' 
#' \item{Bhat}{p x r matrix of the regression coeffcients.}
#'
#' \item{Shat}{p x r matrix of the standard errors for regression coeffcients.}
#' 
#' @export
compute_univariate_sumstats <- function(X, Y, standardize=FALSE, standardize.response=FALSE){
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

###Functions to compute data-driven matrices using flashr
my_init_fn <- function(Y, K = 1) {
  ret <- flashr:::udv_si(Y, K)
  pos_sum <- sum(ret$v[ret$v > 0])
  neg_sum <- -sum(ret$v[ret$v < 0])
  if (neg_sum > pos_sum) {
    return(list(u = -ret$u, d = ret$d, v = -ret$v))
  } else
    return(ret)
}

flash_pipeline <- function(data, ...) {
  ## current state-of-the art
  ## suggested by Jason Willwerscheid
  ## cf: discussion section of
  ## https://willwerscheid.github.io/MASHvFLASH/MASHvFLASHnn2.html
  ebnm_fn <- "ebnm_ash"
  ebnm_param <- list(l = list(mixcompdist = "normal",
                              optmethod = "mixSQP"),
                     f = list(mixcompdist = "+uniform",
                              optmethod = "mixSQP"))
  
  ##
  fl_g <- flashr:::flash_greedy_workhorse(data,
                                          var_type = "constant",
                                          ebnm_fn = ebnm_fn,
                                          ebnm_param = ebnm_param,
                                          init_fn = my_init_fn,
                                          stopping_rule = "factors",
                                          tol = 1e-3,
                                          verbose_output = "odF")
  fl_b <- flashr:::flash_backfit_workhorse(data,
                                           f_init = fl_g,
                                           var_type = "constant",
                                           ebnm_fn = ebnm_fn,
                                           ebnm_param = ebnm_param,
                                           stopping_rule = "factors",
                                           tol = 1e-3,
                                           verbose_output = "odF")
  return(fl_b)
}

cov_flash <- function(data, subset = NULL, non_canonical = FALSE, save_model = NULL) {
  if(is.null(subset)) subset <- 1:mashr:::n_effects(data)
  b.center <- apply(data$Bhat[subset,], 2, function(x) x - mean(x))
  ## Only keep factors with at least two values greater than 1 / sqrt(n)
  find_nonunique_effects <- function(fl) {
    thresh <- 1/sqrt(ncol(fl$fitted_values))
    vals_above_avg <- colSums(fl$ldf$f > thresh)
    nonuniq_effects <- which(vals_above_avg > 1)
    return(fl$ldf$f[, nonuniq_effects, drop = FALSE])
  }
  
  fmodel <- flash_pipeline(b.center)
  if (non_canonical)
    flash_f <- find_nonunique_effects(fmodel)
  else 
    flash_f <- fmodel$ldf$f
  ## row.names(flash_f) <- colnames(b)
  if (!is.null(save_model)) saveRDS(list(model=fmodel, factors=flash_f), save_model)
  if(ncol(flash_f) == 0){
    U.flash <- list("tFLASH" = t(fmodel$fitted_values) %*% fmodel$fitted_values / nrow(fmodel$fitted_values))
  } else{
    U.flash <- c(mashr:::cov_from_factors(t(as.matrix(flash_f)), "FLASH"),
                 list("tFLASH" = t(fmodel$fitted_values) %*% fmodel$fitted_values / nrow(fmodel$fitted_values)))
  }
  return(U.flash)
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


###Compute canonical covariance matrices for mr.mash.
##N.B. deprecated and left only for backwards compatibility
compute_cov_canonical <- function(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1), grid, zeromat=TRUE){
  S <- create_cov_canonical(r, singletons, hetgrid)
  U <- list()
  nms <- vector("character")
  t <- 0
  for(i in 1:length(S)){
    for(j in 1:length(grid)){
      t <- t+1
      U[[t]] <- S[[i]]*grid[j]
      nms[t] <- paste0(names(S)[i], "_grid", j) 
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

###Compute univariate summary statistics
##N.B. deprecated and left only for backwards compatibility
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

