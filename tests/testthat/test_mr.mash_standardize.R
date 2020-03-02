context("Test mr.mash with standardize=T vs old mr.mash.scaled.X")

test_that("mr.mash with standardize=T and old mr.mash.scaled.X return the same results", {
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
  
  ###Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance
  Y <- sim_mvr(X, B, V)
  
  ###Specify the mixture weights and covariance matrices for the mixture-of-normals prior
  grid <- seq(1, 5)
  S0mix <- compute_cov_canonical(ncol(Y), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid, zeromat=TRUE)
  
  w0    <- rep(1/(length(S0mix)), length(S0mix))
  
  ###Estimate residual covariance
  V_est <- cov(Y)
  
  ###Fit with current implementation
  fit_current <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE, update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE, verbose=FALSE, update_V=FALSE)
  
  ###Load old fit
  fit_old <- readRDS("mr.mash.scaled.X_fit.rds")
  
  ###Tests
  expect_equal(fit_current$mu1, fit_old$mu1, tolerance = 1e-10, scale = 1)
  expect_equal(fit_current$S1, fit_old$S1, tolerance = 1e-10, scale = 1)
  expect_equal(fit_current$w1, fit_old$w1, tolerance = 1e-10, scale = 1)
  expect_equal(fit_current$ELBO, fit_old$ELBO, tolerance = 1e-10, scale = 1)
})
