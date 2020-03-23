context("Test current vs old mr.mash implementations")

test_that("Current and old mr.mash return the same results",{
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
  grid  <- seq(1, 5)
  S0mix <- compute_cov_canonical(ncol(Y), singletons=TRUE,
                                 hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid,
                                 zeromat=TRUE)
  w0    <- rep(1/(length(S0mix)), length(S0mix))
  
  ###Estimate residual covariance
  V_est <- cov(Y)
  
  ###Fit with current implementation
  capture.output(
    fit_current <- mr.mash(X, Y, S0mix, w0, V_est, tol=1e-8, update_w0=TRUE,
                           update_w0_method="EM", compute_ELBO=TRUE, 
                           standardize=FALSE, verbose=FALSE, update_V=FALSE,
                           version="R"))
  capture.output(
    fit_current_rcpp <- mr.mash(X, Y, S0mix, w0, V_est, tol=1e-8,
                                update_w0=TRUE, update_w0_method="EM",
                                compute_ELBO=TRUE, standardize=FALSE,
                                verbose=FALSE, update_V=FALSE, version="Rcpp"))
  
  ###Load old fit
  fit_old <- readRDS("mr_mash_slow_implement_fit.rds")
  
  ###Tests
  expect_equal(fit_current$mu1, fit_old$mu1, tolerance = 1e-6, scale = 1)
  expect_equal(fit_current$S1, fit_old$S1, tolerance = 1e-6, scale = 1)
  expect_equal(fit_current$w1, fit_old$w1, tolerance = 1e-6, scale = 1)
  expect_equal(fit_current$ELBO, fit_old$ELBO, tolerance = 9.99e-5, scale = 1)
  
  expect_equal(fit_current_rcpp$mu1, fit_old$mu1, tolerance = 1e-6, scale = 1)
  expect_equal(fit_current_rcpp$S1, fit_old$S1, tolerance = 1e-6, scale = 1)
  expect_equal(fit_current_rcpp$w1, fit_old$w1, tolerance = 1e-6, scale = 1)
  expect_equal(fit_current_rcpp$ELBO, fit_old$ELBO, tolerance = 9.99e-5,
               scale = 1)
  
})
