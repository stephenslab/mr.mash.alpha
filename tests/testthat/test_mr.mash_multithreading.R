context("mr.mash with 1 or 2 thread(s) return the same results")

test_that("mr.mash with 1 or 2 thread(s) return the same results", {
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
  
  ###Compute quantities needed for mr.mash.rss
  out <- compute_univariate_sumstats(X=X, Y=Y, standardize=FALSE, standardize.response=FALSE, mc.cores=1)
  R <- cor(X)
  X_colMeans <- colMeans(X)
  Y_colMeans <- colMeans(Y)
  
  ###Fit with current implementation (1 thread)
  capture.output(
    fit_1 <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                   update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE,
                   verbose=FALSE, update_V=TRUE, nthreads=1))
  fit_1$progress <- fit_1$progress[, -2] ##This line is needed to remove the timing column -->  
                                         ##hopefully faster when using multiple threads
  capture.output(
    fit_1_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE,
                          verbose=FALSE, update_V=TRUE, nthreads=1))
  fit_1_miss$progress <- fit_1_miss$progress[, -2] ##This line is needed to remove the timing column -->  
  ##hopefully faster when using multiple threads
  
  fit_1_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                           w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, standardize=TRUE,
                           verbose=FALSE, update_V=TRUE, X_colmeans=X_colMeans, Y_colmeans=Y_colMeans, 
                           nthreads=1)
  fit_1_rss$progress <- fit_1_rss$progress[, -2] ##This line is needed to remove the timing column -->  
  ##hopefully faster when using multiple threads

  
  
  ###Fit with current implementation (2 threads)
  capture.output(
    fit_2 <- mr.mash(X, Y, S0mix, w0, V_est, update_w0=TRUE,
                     update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE,
                     verbose=FALSE, update_V=TRUE, nthreads=2))
  fit_2$progress <- fit_2$progress[, -2]
  
  capture.output(
    fit_2_miss <- mr.mash(X, Y_miss, S0mix, w0, V_est, update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE,
                          verbose=FALSE, update_V=TRUE, nthreads=2))
  fit_2_miss$progress <- fit_2_miss$progress[, -2]
  
  fit_2_rss <- mr.mash.rss(Bhat=out$Bhat, Shat=out$Shat, covY=V_est, R=R, n=n, S0=S0mix, 
                           w0=w0, V=V_est, update_w0=TRUE, compute_ELBO=TRUE, standardize=TRUE,
                           verbose=FALSE, update_V=TRUE, X_colmeans=X_colMeans, Y_colmeans=Y_colMeans, 
                           nthreads=2)
  fit_2_rss$progress <- fit_2_rss$progress[, -2] ##This line is needed to remove the timing column -->  
  ##hopefully faster when using multiple threads
  
  
  
  ###Test
  expect_equal(fit_1, fit_2, tolerance=1e-10, scale=1)
  expect_equal(fit_1_miss, fit_2_miss, tolerance=1e-10, scale=1)
})