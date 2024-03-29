context("mr.mash R and Rcpp versions return same result")

test_that("mr.mash R version and Rcpp version return the same results", {
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
  
  ###Assign some missing values in Y
  Y_miss <- Y
  Y_miss[1:5, 1] <- NA
  Y_miss[6:10, 2] <- NA
  
  ###Specify the mixture weights and covariance matrices for the
  ###mixture-of-normals prior
  grid <- seq(1, 5)
  S0mix <- compute_cov_canonical(ncol(Y), singletons=TRUE,
                                 hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid,
                                 zeromat=TRUE)
  
  w0    <- rep(1/(length(S0mix)), length(S0mix))
  
  ###Estimate residual covariance
  V_est <- cov(Y)
  
  ###Fit with current implementation (R version)
  capture.output(
    fit <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                   update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE,
                   verbose=FALSE, update_V=FALSE, version="R"))
  fit$progress <- fit$progress[, -2] ##This line is needed to remove the timing column --> clearly different between R and Rcpp
  
  capture.output(
    fit_scaled <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, 
                          standardize=TRUE, verbose=FALSE, update_V=FALSE,
                          version="R"))
  fit_scaled$progress <- fit_scaled$progress[, -2]
  
  capture.output(
    fit_V <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                     update_w0_method="EM", compute_ELBO=TRUE, 
                     standardize=FALSE, verbose=FALSE, update_V=TRUE,
                     version="R"))
  fit_V$progress <- fit_V$progress[, -2]
  
  capture.output(
    fit_scaled_V <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                            update_w0_method="EM", compute_ELBO=TRUE, 
                            standardize=TRUE, verbose=FALSE, update_V=TRUE,
                            version="R"))
  fit_scaled_V$progress <- fit_scaled_V$progress[, -2]
  
  capture.output(
    fit_scaled_V_declogBF <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                            update_w0_method="EM", compute_ELBO=TRUE, 
                            standardize=TRUE, verbose=FALSE, update_V=TRUE,
                            version="R", ca_update_order="decreasing_logBF"))
  fit_scaled_V_declogBF$progress <- fit_scaled_V_declogBF$progress[, -2]

  capture.output(
    fit_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                   update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE,
                   verbose=FALSE, update_V=FALSE, version="R"))
  fit_miss$progress <- fit_miss$progress[, -2]
  
  capture.output(
    fit_scaled_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, 
                          standardize=TRUE, verbose=FALSE, update_V=FALSE,
                          version="R"))
  fit_scaled_miss$progress <- fit_scaled_miss$progress[, -2]
  
  capture.output(
    fit_V_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                     update_w0_method="EM", compute_ELBO=TRUE, 
                     standardize=FALSE, verbose=FALSE, update_V=TRUE,
                     version="R"))
  fit_V_miss$progress <- fit_V_miss$progress[, -2]
  
  capture.output(
    fit_scaled_V_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                            update_w0_method="EM", compute_ELBO=TRUE, 
                            standardize=TRUE, verbose=FALSE, update_V=TRUE,
                            version="R"))
  fit_scaled_V_miss$progress <- fit_scaled_V_miss$progress[, -2]
  

    ###Fit with current implementation (Rcpp version)
  capture.output(
    fit_rcpp <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                        update_w0_method="EM", compute_ELBO=TRUE, 
                        standardize=FALSE, verbose=FALSE, update_V=FALSE,
                        version="Rcpp"))
  fit_rcpp$progress <- fit_rcpp$progress[, -2]
  
  capture.output(
    fit_scaled_rcpp <- mr.mash(X, Y, S0mix, w0, V_est,
                               update_w0=TRUE, update_w0_method="EM",
                               compute_ELBO=TRUE, standardize=TRUE,
                               verbose=FALSE, update_V=FALSE, version="Rcpp"))
  fit_scaled_rcpp$progress <- fit_scaled_rcpp$progress[, -2]
  
  capture.output(
    fit_V_rcpp <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, 
                          standardize=FALSE, verbose=FALSE, update_V=TRUE,
                          version="Rcpp"))
  fit_V_rcpp$progress <- fit_V_rcpp$progress[, -2]
  
  capture.output(
    fit_scaled_V_rcpp <- mr.mash(X, Y, S0mix, w0, V_est,
                                 update_w0=TRUE, update_w0_method="EM",
                                 compute_ELBO=TRUE, standardize=TRUE,
                                 verbose=FALSE, update_V=TRUE, version="Rcpp"))
  fit_scaled_V_rcpp$progress <- fit_scaled_V_rcpp$progress[, -2]
  
  capture.output(
    fit_scaled_V_declogBF_rcpp <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                                     update_w0_method="EM", compute_ELBO=TRUE, 
                                     standardize=TRUE, verbose=FALSE, update_V=TRUE,
                                     version="Rcpp", ca_update_order="decreasing_logBF"))
  fit_scaled_V_declogBF_rcpp$progress <- fit_scaled_V_declogBF_rcpp$progress[, -2]
  
  capture.output(
    fit_miss_rcpp <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                        update_w0_method="EM", compute_ELBO=TRUE, 
                        standardize=FALSE, verbose=FALSE, update_V=FALSE,
                        version="Rcpp"))
  fit_miss_rcpp$progress <- fit_miss_rcpp$progress[, -2]
  
  capture.output(
    fit_scaled_miss_rcpp <- mr.mash(X, Y_miss, S0mix, w0, V_est,
                               update_w0=TRUE, update_w0_method="EM",
                               compute_ELBO=TRUE, standardize=TRUE,
                               verbose=FALSE, update_V=FALSE, version="Rcpp"))
  fit_scaled_miss_rcpp$progress <- fit_scaled_miss_rcpp$progress[, -2]
  
  capture.output(
    fit_V_miss_rcpp <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, 
                          standardize=FALSE, verbose=FALSE, update_V=TRUE,
                          version="Rcpp"))
  fit_V_miss_rcpp$progress <- fit_V_miss_rcpp$progress[, -2]
  
  capture.output(
    fit_scaled_V_miss_rcpp <- mr.mash(X, Y_miss, S0mix, w0, V_est,
                                 update_w0=TRUE, update_w0_method="EM",
                                 compute_ELBO=TRUE, standardize=TRUE,
                                 verbose=FALSE, update_V=TRUE, version="Rcpp"))
  fit_scaled_V_miss_rcpp$progress <- fit_scaled_V_miss_rcpp$progress[, -2]
  
  
  
  
  ###Tests
  expect_equal(fit, fit_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled, fit_scaled_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_V, fit_V_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled_V, fit_scaled_V_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled_V_declogBF, fit_scaled_V_declogBF_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_miss, fit_miss_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled_miss, fit_scaled_miss_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_V_miss, fit_V_miss_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled_V_miss, fit_scaled_V_miss_rcpp, tolerance=1e-10, scale=1)
})
