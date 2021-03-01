context("mr.mash and mr_mash_simple return same result")

test_that("mr.mash and mr_mash_simple return the same results", {
  ###Set seed
  set.seed(123)
  
  ###Simulate X and Y
  n  <- 100
  p <- 10
  r <- 2
  
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
  
  ###Initial value for the regression coefficients
  B0  <- matrix(0,p,r)
  
  ###Estimate residual covariance
  V_est <- cov(Y)
  
  ###Fit with mr.mash (R version)
  capture.output(
    fit <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                  update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE,
                  verbose=FALSE, update_V=FALSE, version="R", mu1_init=B0))
  
  capture.output(
    fit_V <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                    update_w0_method="EM", compute_ELBO=TRUE, 
                    standardize=FALSE, verbose=FALSE, update_V=TRUE,
                    version="R", mu1_init=B0))
  
  capture.output(
    fit_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                   update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE,
                   verbose=FALSE, update_V=FALSE, version="R", mu1_init=B0))

  capture.output(
    fit_V_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                     update_w0_method="EM", compute_ELBO=TRUE, 
                     standardize=FALSE, verbose=FALSE, update_V=TRUE,
                     version="R", mu1_init=B0))

  
  ###Fit with mr_mash_simple
  capture.output(
    fit_simple <- mr_mash_simple(X, Y, V_est, lapply(S0mix, makePD, e=1e-8), w0, update_w0=TRUE,
                                 verbose=FALSE, update_V=FALSE, B=B0))
  
  capture.output(
    fit_V_simple <- mr_mash_simple(X, Y, V_est, lapply(S0mix, makePD, e=1e-8), w0, update_w0=TRUE,
                                   verbose=FALSE, update_V=TRUE, B=B0))
  
  capture.output(
    fit_miss_simple <- mr_mash_simple(X, Y_miss, V_est, lapply(S0mix, makePD, e=1e-8), w0, update_w0=TRUE,
                                      verbose=FALSE, update_V=FALSE, B=B0))
  
  capture.output(
    fit_V_miss_simple <- mr_mash_simple(X, Y_miss, V_est, lapply(S0mix, makePD, e=1e-8), w0, update_w0=TRUE,
                                        verbose=FALSE, update_V=TRUE, B=B0))
  

  
  ###Tests
  expect_equivalent(fit$intercept, fit_simple$intercept, tolerance=1e-10, scale=1)
  expect_equivalent(fit$mu1, fit_simple$B, tolerance=1e-10, scale=1)
  expect_equivalent(fit$w0, fit_simple$w0, tolerance=1e-10, scale=1)
  expect_equivalent(fit$progress$ELBO, fit_simple$ELBO, tolerance=1e-10, scale=1)
  
  expect_equivalent(fit_V$intercept, fit_V_simple$intercept, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V$mu1, fit_V_simple$B, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V$w0, fit_V_simple$w0, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V$progress$ELBO, fit_V_simple$ELBO, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V$V, fit_V_simple$V, tolerance=1e-10, scale=1)
  
  
  expect_equivalent(fit_miss$intercept, fit_miss_simple$intercept, tolerance=1e-10, scale=1)
  expect_equivalent(fit_miss$mu1, fit_miss_simple$B, tolerance=1e-10, scale=1)
  expect_equivalent(fit_miss$w0, fit_miss_simple$w0, tolerance=1e-10, scale=1)
  expect_equivalent(fit_miss$progress$ELBO, fit_miss_simple$ELBO, tolerance=1e-10, scale=1)
  expect_equivalent(fit_miss$Y, fit_miss_simple$Y, tolerance=1e-10, scale=1)
  
  expect_equivalent(fit_V_miss$intercept, fit_V_miss_simple$intercept, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V_miss$mu1, fit_V_miss_simple$B, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V_miss$w0, fit_V_miss_simple$w0, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V_miss$progress$ELBO, fit_V_miss_simple$ELBO, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V_miss$Y, fit_V_miss_simple$Y, tolerance=1e-10, scale=1)
  expect_equivalent(fit_V_miss$V, fit_V_miss_simple$V, tolerance=1e-10, scale=1)
})
