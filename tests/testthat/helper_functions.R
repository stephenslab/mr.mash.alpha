# Used to check that x is a vector in which x[i+1] >= x[i] for all i.
expect_nondecreasing <- function(x, tolerance, scale){
  expect_equal(diff(x) >= 0,rep(TRUE,length(x) - 1), tolerance = tolerance, scale = scale)
}
  
