context("Test Bayesian multivariate regression with transformed x")

test_that("bayes_mvr_mix and bayes_mvr_mix_transformed_X return the same results", {
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
                 rep(0, (p-2)*2)), byrow=T, ncol=2)
  
  ###Simulate X
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  X <- scale(X, center=T, scale=F)
  
  ###Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance
  Y <- sim_mvr(X, B, V)
  
  ###Specify the mixture weights and covariance matrices for the mixture-of-normals prior
  grid <- seq(1, 5)
  S0mix <- compute_cov_canonical(ncol(Y), singletons=T, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid, zeromat=T)
  
  w0    <- rep(1/(length(S0mix)), length(S0mix))
  
  ###Fit BMR with untransformed X
  fit_mix <- bayes_mvr_mix(X[, 1], Y, V, w0, S0mix)
  
  ###Fit BMR with transformed X
  comps1 <- precompute_quants_transformed_X(X, V, S0mix)
  fit_mix_transformed <- bayes_mvr_mix_transformed_X(X[, 1], Y, V, w0, S0mix, comps1$xtx[1], comps1$V_chol, comps1$U0, comps1$d, comps1$Q)
  
  ###Tests
  expect_equal(fit_mix$mu1, fit_mix_transformed$mu1, tolerance = 1e-10, scale = 1)
  expect_equal(fit_mix$S1, fit_mix_transformed$S1, tolerance = 1e-10, scale = 1)
  expect_equal(fit_mix$w1, fit_mix_transformed$w1, tolerance = 1e-10, scale = 1)
  expect_equal(fit_mix$logbf, fit_mix_transformed$logbf, tolerance = 1e-10, scale = 1)
})
