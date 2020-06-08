context("Test scale_fast")

test_that("scale and scale_fast return the same result",{
  set.seed(123)
  
  ###Simulate data
  n <- 1000
  p <- 1000
  M <- matrix(rnorm(n*p), ncol=p, nrow=n)
  colnames(M) <- paste0("X", seq(1, p))
  rownames(M) <- paste0("N", seq(1, n))
  
  ###Assign missing values
  M[sample(size=100, x=seq(n*p))] <- NA
  
  ###Scale M
  Ms_fast <- scale_fast(M, scale=TRUE)
  Ms_fast2 <- scale_fast2(M, scale=TRUE)
  Ms <- scale(M, scale=TRUE)
  Mc_fast <- scale_fast(M, scale=FALSE)
  Mc_fast2 <- scale_fast2(M, scale=FALSE)
  Mc <- scale(M, scale=FALSE)
  
  Ms_means <- attr(Ms,"scaled:center")
  Ms_sds <- attr(Ms,"scaled:scale")
  
  attr(Ms,"scaled:center") <- NULL
  attr(Ms,"scaled:scale") <- NULL
  
  Mc_means <- attr(Mc,"scaled:center")
  
  attr(Mc,"scaled:center") <- NULL
  
  ###Tests
  expect_equal(Ms_fast$M, Ms, tolerance = 1e-10, scale = 1)
  expect_equal(Ms_fast$means, Ms_means, tolerance = 1e-10, scale = 1)
  expect_equivalent(Ms_fast$sds, Ms_sds, tolerance = 1e-10, scale = 1)
  
  expect_equal(Ms_fast2$M, Ms, tolerance = 1e-10, scale = 1)
  expect_equal(Ms_fast2$means, Ms_means, tolerance = 1e-10, scale = 1)
  expect_equivalent(Ms_fast2$sds, Ms_sds, tolerance = 1e-10, scale = 1)
  
  expect_equal(Mc_fast$M, Mc, tolerance = 1e-10, scale = 1)
  expect_equal(Mc_fast$means, Mc_means, tolerance = 1e-10, scale = 1)
  
  expect_equal(Mc_fast2$M, Mc, tolerance = 1e-10, scale = 1)
  expect_equal(Mc_fast2$means, Mc_means, tolerance = 1e-10, scale = 1)
})