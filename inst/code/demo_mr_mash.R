# An illustration of the mr_mash_simple implementation applied to a
# small, simulated data set.
suppressMessages(library(MBSP))
library(mvtnorm)

# SCRIPT PARAMETERS
# -----------------
# Number of samples (n) and number of predictors (p).
n <- 500
p <- 20

# Residual covariance matrix.
V <- rbind(c(1.0,0.2),
           c(0.2,0.4))
r <- nrow(V)
rownames(V) <- paste0("r",1:r)
colnames(V) <- paste0("r",1:r)

# True effects used to simulate the data.
intercept <- c(-1,+1)
B <- rbind(c(-2.0, -1.5),
           c( 1.0,  1.0),
           matrix(0,p - 2,r))
rownames(B) <- paste0("x",1:p)
colnames(B) <- paste0("d",1:r)

# Covariances in the mixture-of-normals prior on the regression
# coefficients.
S0 <- list(k1 = rbind(c(3,0),
                      c(0,3)),
           k2 = rbind(c(4,2),
                      c(2,4)),
           k3 = rbind(c(6,3.5),
                      c(3.5,4)),
           k4 = rbind(c(5,0),
                      c(0,0)))
S0 <- lapply(S0,function (x) {
  rownames(x) <- paste0("r",1:r)
  colnames(x) <- paste0("r",1:r)
  return(x)
})

# The mixture weights in the mixture-of-normals prior on the
# regression coefficients.
w0 <- c(0.1,0.6,0.2,0.1)
k  <- length(w0)
names(w0) <- paste0("k",1:k)

# SIMULATE DATA
# -------------
set.seed(1)
X <- matrix(rnorm(n*p),n,p)
rownames(X) <- paste0("i",1:n)
colnames(X) <- paste0("x",1:p)

# Simulate Y ~ MN(X*B,I,V). Note that matrix.normal from the MBSP
# package appears to be much faster than rmatrixnorm from the
# MixMatrix package.
Y <- matrix.normal(X %*% B + matrix(intercept,n,r,byrow = TRUE),diag(n),V)

# FIT MR-MASH MODEL
# -----------------
# Run 50 co-ordinate ascent updates.
mu1 <- matrix(0,p,r)
fit <- mr.mash(X,Y,V,S0,w0,mu_init = mu1,max_iter = 50,
               standardize = FALSE,version = "Rcpp")

# Compare the posterior mean estimates of the regression coefficients
# against the coefficients used to simulate the data.
plot(B,fit$mu1,pch = 20,xlab = "true",ylab = "estimated")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

# Compare the fitted responses against the ground-truth responses.
Yest <- predict(fit,X)
plot(Y,Yest,pch = 20)
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")
