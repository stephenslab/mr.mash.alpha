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
    fit <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE,
                   update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE,
                   verbose=FALSE, update_V=FALSE, version="R"))
  capture.output(
    fit_scaled <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, 
                          standardize=TRUE, verbose=FALSE, update_V=FALSE,
                          version="R"))
  capture.output(
    fit_V <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE,
                     update_w0_method="EM", compute_ELBO=TRUE, 
                     standardize=FALSE, verbose=FALSE, update_V=TRUE,
                     version="R"))
  capture.output(
    fit_scaled_V <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE,
                            update_w0_method="EM", compute_ELBO=TRUE, 
                            standardize=TRUE, verbose=FALSE, update_V=TRUE,
                            version="R"))
  
  capture.output(
    fit_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7, update_w0=TRUE,
                          update_w0_method="mixsqp", compute_ELBO=TRUE, 
                          standardize=FALSE, verbose=FALSE, update_V=FALSE,
                          version="R"))
  capture.output(
    fit_scaled_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7,
                                 update_w0=TRUE, update_w0_method="mixsqp",
                                 compute_ELBO=TRUE, standardize=TRUE,
                                 verbose=FALSE, update_V=FALSE, version="R"))
  capture.output(
    fit_V_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7, update_w0=TRUE,
                            update_w0_method="mixsqp", compute_ELBO=TRUE, 
                            standardize=FALSE, verbose=FALSE, update_V=TRUE,
                            version="R"))
  capture.output(
    fit_scaled_V_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7,
                                   update_w0=TRUE, update_w0_method="mixsqp",
                                   compute_ELBO=TRUE, standardize=TRUE,
                                   verbose=FALSE, update_V=TRUE, version="R"))
  
  ###Fit with current implementation (Rcpp version)
  capture.output(
    fit_rcpp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE,
                        update_w0_method="EM", compute_ELBO=TRUE, 
                        standardize=FALSE, verbose=FALSE, update_V=FALSE,
                        version="Rcpp"))
  capture.output(
    fit_scaled_rcpp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8,
                               update_w0=TRUE, update_w0_method="EM",
                               compute_ELBO=TRUE, standardize=TRUE,
                               verbose=FALSE, update_V=FALSE, version="Rcpp"))
  capture.output(
    fit_V_rcpp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, 
                          standardize=FALSE, verbose=FALSE, update_V=TRUE,
                          version="Rcpp"))
  capture.output(
    fit_scaled_V_rcpp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8,
                                 update_w0=TRUE, update_w0_method="EM",
                                 compute_ELBO=TRUE, standardize=TRUE,
                                 verbose=FALSE, update_V=TRUE, version="Rcpp"))
  
  capture.output(
    fit_rcpp_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7,
                               update_w0=TRUE, update_w0_method="mixsqp",
                               compute_ELBO=TRUE, standardize=FALSE,
                               verbose=FALSE, update_V=FALSE, version="Rcpp"))
  capture.output(
    fit_scaled_rcpp_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7,
                                      update_w0=TRUE,
                                      update_w0_method="mixsqp",
                                      compute_ELBO=TRUE, standardize=TRUE,
                                      verbose=FALSE, update_V=FALSE,
                                      version="Rcpp"))
  capture.output(
    fit_V_rcpp_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7,
                                 update_w0=TRUE, update_w0_method="mixsqp",
                                 compute_ELBO=TRUE, standardize=FALSE,
                                 verbose=FALSE, update_V=TRUE, version="Rcpp"))
  capture.output(
    fit_scaled_V_rcpp_mixsqp <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-7,
                                        update_w0=TRUE,
                                        update_w0_method="mixsqp",
                                        compute_ELBO=TRUE, standardize=TRUE,
                                        verbose=FALSE, update_V=TRUE,
                                        version="Rcpp"))
  
  ###Tests
  expect_equal(fit, fit_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled, fit_scaled_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_V, fit_V_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled_V, fit_scaled_V_rcpp, tolerance=1e-10, scale=1)
  expect_equal(fit_mixsqp, fit_rcpp_mixsqp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled_mixsqp, fit_scaled_rcpp_mixsqp, tolerance=1e-10,
               scale=1)
  expect_equal(fit_V_mixsqp, fit_V_rcpp_mixsqp, tolerance=1e-10, scale=1)
  expect_equal(fit_scaled_V_mixsqp, fit_scaled_V_rcpp_mixsqp, tolerance=1e-10,
               scale=1)
})
