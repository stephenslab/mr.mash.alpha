context("Test mr.mash outputted fitted values vs predicted values with the same X")

test_that("mr.mash outputted fitted values vs predicted values with the same X are equal", {
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
  
  ###Estimate residual covariance
  V_est <- cov(Y)
  
  ###Fit the model
  fit <- mr.mash(Y, X, V_est, S0mix, w0, tol=1e-8, update_w0=T, compute_ELBO=T, standardize=T, verbose=F)
  
  ###Predict values with tha same X 
  Yhat <- predict(fit, X)
  
  ###Tests
  expect_equal(fit$fitted, Yhat, tolerance = 1e-10, scale = 1)
})
