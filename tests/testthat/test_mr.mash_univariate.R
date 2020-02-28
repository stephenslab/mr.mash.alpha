context("Test mr.mash with univariate Y")

test_that("mr.mash and varbvsmix return the same results with univariate Y", {
  ###Set options
  options(stringsAsFactors = FALSE)
  
  ###Set seed
  set.seed(123)
  
  ###Simulate X and Y
  ##Set parameters
  n  <- 100
  p <- 10
  
  ##Set residual covariance
  V  <- matrix(1.4, ncol=1, nrow=1)
  
  ##Set true effects
  B  <- matrix(c(-2, 5, rep(0, (p-2))), byrow=TRUE, ncol=1)
  
  ##Simulate X
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  X <- scale(X, center=TRUE, scale=FALSE)
  
  ##Simulate Y
  Y <- X%*%B + rnorm(n=n, sd=sqrt(as.numeric(V)))
  
  ###Specify the mixture weights and covariance matrices for the mixture-of-normals prior.
  grid <- seq(0, 5)
  
  S0mix <- list()
  for(i in 1:length(grid)){
    S0mix[[i]] <- matrix(grid[i], ncol=1, nrow=1)
  }
  
  w0    <- rep(1/(length(S0mix)), length(S0mix))
  
  ###Estimate residual covariance
  V_est <- cov(Y)
  
  ###Fit mr.mash
  fit_mr.mash <- mr.mash(X, Y, V_est, S0mix, w0, tol=1e-8, update_w0=TRUE, update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE)
  
  ###Fit varbvsmix
  fit_varbvsmix <- varbvs::varbvsmix(X,NULL,Y,sigma = as.numeric(V_est),w = w0,
                                     sa = grid/as.numeric(V_est),
                                     update.sigma = FALSE,update.w = TRUE,
                                     drop.threshold = 0,tol = 1e-8,
                                     verbose = FALSE)
  mu1_varbvsmix <- rowSums(with(fit_varbvsmix,alpha * mu))
  betavarmix <- function (p, mu, s){
    rowSums(p*(s + mu^2)) - rowSums(p*mu)^2  
  }
  S1_varbvsmix <- betavarmix(fit_varbvsmix$alpha, fit_varbvsmix$mu, fit_varbvsmix$s)
  
  ###Tests
  expect_equivalent(drop(fit_mr.mash$mu1), mu1_varbvsmix, tolerance = 1e-6, scale = 1)
  expect_equivalent(drop(fit_mr.mash$S1), S1_varbvsmix, tolerance = 1e-6, scale = 1)
  expect_equivalent(fit_mr.mash$w1, fit_varbvsmix$alpha, tolerance = 1e-6, scale = 1)
  expect_equivalent(fit_mr.mash$ELBO, tail(fit_varbvsmix$logZ, 1), tolerance = 1e-6, scale = 1)
})
