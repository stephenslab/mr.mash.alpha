context("mr.mash and mr.mash.rss versions return same result")

test_that("mr.mash and mr.mash.rss return the same results", {
  ###Set seed
  set.seed(123)
  
  ###Simulate X and Y
  n  <- 100
  p <- 10
  
  ###Set residual covariance
  V  <- rbind(c(1.0,0.2),
              c(0.2,0.4))
  
  ###Set true effects
  B  <- matrix(c(-2, -2,
                 5, 5,
                 rep(0, (p-2)*2)), byrow=TRUE, ncol=2)
  
  ###Simulate X
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  X <- scale(X, center=TRUE, scale=FALSE)
  
  ###Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity
  ###matrix and V is the residual covariance
  Y <- sim_mvr(X, B, V)
  
  ###Specify the mixture weights and covariance matrices for the
  ###mixture-of-normals prior
  grid <- seq(1, 5)
  S0mix <- compute_cov_canonical(ncol(Y), singletons=TRUE,
                                 hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid,
                                 zeromat=TRUE)
  
  w0    <- rep(1/(length(S0mix)), length(S0mix))
  
  ###Estimate residual covariance
  V_est <- cov(Y)
  
  ###Fit mr.mash
  fit <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                   update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE,
                   verbose=FALSE, update_V=FALSE)
  fit$progress <- fit$fitted <- fit$pve <- fit$G <- NULL
  
  fit_scaled <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, 
                          standardize=TRUE, verbose=FALSE, update_V=FALSE)
  fit_scaled$progress <- fit_scaled$fitted <- fit_scaled$pve <- fit_scaled$G <- NULL
  
  fit_V <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                     update_w0_method="EM", compute_ELBO=TRUE, 
                     standardize=FALSE, verbose=FALSE, update_V=TRUE)
  fit_V$progress <- fit_V$fitted <- fit_V$pve <- fit_V$G <- NULL
  
  fit_V_diag <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                   update_w0_method="EM", compute_ELBO=TRUE, 
                   standardize=FALSE, verbose=FALSE, update_V=TRUE,
                   update_V_method="diagonal")
  fit_V_diag$progress <- fit_V_diag$fitted <- fit_V_diag$pve <- fit_V_diag$G <- NULL
  
  fit_scaled_V <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                            update_w0_method="EM", compute_ELBO=TRUE, 
                            standardize=TRUE, verbose=FALSE, update_V=TRUE)
  fit_scaled_V$progress <- fit_scaled_V$fitted <- fit_scaled_V$pve <- fit_scaled_V$G <- NULL
  
  fit_scaled_V_declogBF <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                            update_w0_method="EM", compute_ELBO=TRUE, 
                            standardize=TRUE, verbose=FALSE, update_V=TRUE,
                            ca_update_order="decreasing_logBF")
  fit_scaled_V_declogBF$progress <- fit_scaled_V_declogBF$fitted <- fit_scaled_V_declogBF$pve <- fit_scaled_V_declogBF$G <- NULL


  ###Fit mr.mash.rss
  out <- compute_univariate_sumstats(X=X, Y=Y, standardize=FALSE, standardize.response=FALSE, mc.cores=1)
  R <- cor(X)
  X_colMeans <- colMeans(X)
  Y_colMeans <- colMeans(Y)
  
  fit_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                          w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, standardize=FALSE,
                          verbose=FALSE, update_V=FALSE, X_colmeans=X_colMeans, Y_colmeans=Y_colMeans)
  fit_rss$progress <- NULL
  
  fit_scaled_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                                w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, 
                                standardize=TRUE, verbose=FALSE, update_V=FALSE,
                                X_colmeans=X_colMeans, Y_colmeans=Y_colMeans)
  fit_scaled_rss$progress <- NULL
  
  fit_V_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                          w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, 
                          standardize=FALSE, verbose=FALSE, update_V=TRUE,
                          X_colmeans=X_colMeans, Y_colmeans=Y_colMeans)
  fit_V_rss$progress <- NULL
  
  fit_V_diag_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                           w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, 
                           standardize=FALSE, verbose=FALSE, update_V=TRUE,
                           X_colmeans=X_colMeans, Y_colmeans=Y_colMeans,
                           update_V_method="diagonal")
  fit_V_diag_rss$progress <- NULL
  
  fit_scaled_V_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                                  w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, 
                                  standardize=TRUE, verbose=FALSE, update_V=TRUE,
                                  X_colmeans=X_colMeans, Y_colmeans=Y_colMeans)
  fit_scaled_V_rss$progress <- NULL
  
  
  fit_scaled_V_declogBF_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                                          w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, 
                                          standardize=TRUE, verbose=FALSE, update_V=TRUE,
                                          ca_update_order="decreasing_logBF",
                                          X_colmeans=X_colMeans, Y_colmeans=Y_colMeans)
  fit_scaled_V_declogBF_rss$progress <- NULL
  
  
  
  
  ###Tests
  expect_equal(unclass(fit), unclass(fit_rss), tolerance=1e-10, scale=1)
  expect_equal(unclass(fit_scaled), unclass(fit_scaled_rss), tolerance=1e-10, scale=1)
  expect_equal(unclass(fit_V), unclass(fit_V_rss), tolerance=1e-10, scale=1)
  expect_equal(unclass(fit_V_diag), unclass(fit_V_diag_rss), tolerance=1e-10, scale=1)
  expect_equal(unclass(fit_scaled_V), unclass(fit_scaled_V_rss), tolerance=1e-10, scale=1)
  expect_equal(unclass(fit_scaled_V_declogBF), unclass(fit_scaled_V_declogBF_rss), tolerance=1e-10, scale=1)
})
