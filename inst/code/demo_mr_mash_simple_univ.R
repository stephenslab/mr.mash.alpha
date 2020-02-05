# A test that mr_mash_simple also works for univariate linear
# regression (i.e., r = 1).
library(mvtnorm)
library(varbvs)
source("../../R/misc.R")
source("../../R/mr_mash_simple.R")

# SCRIPT PARAMETERS
# -----------------
# Number of samples (n) and number of predictors (p).
n <- 500
p <- 20

# True effects used to simulate the data.
b <- c(-2,1,rep(0,p - 2))

# Variances in the mixture-of-normals prior on the regression
# coefficients.
s0 <- list(k1 = 1e-10,k2 = 4,k3 = 6,k4 = 5)

# The mixture weights in the mixture-of-normals prior on the
# regression coefficients.
w0 <- c(0.1,0.6,0.2,0.1)
k  <- length(w0)

# SIMULATE DATA
# -------------
set.seed(1)
X <- matrix(rnorm(n*p),n,p)
X <- scale(X,scale = FALSE)

# Simulate Y ~ MN(X*B,I,V). Note that matrix.normal from the MBSP
# package appears to be much faster than rmatrixnorm from the
# MixMatrix package.
y <- drop(X %*% b + rnorm(n))
y <- y - mean(y)

# FIT MR-MASH MODEL
# -----------------
# Run 20 co-ordinate ascent updates.
b0  <- rep(0,p)
fit <- mr_mash_simple(X,y,1,s0,w0,b0,20)

# Compare the posterior mean estimates of the regression coefficients
# against the coefficients used to simulate the data.
plot(b,fit$B,pch = 20,xlab = "true",ylab = "estimated")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

# FIR VARBVSMIX MODEL
# -------------------
s0    <- unlist(s0)
s0[1] <- 0
out <- varbvsmix(X,NULL,y,s0,1,w0,matrix(0,p,k),matrix(0,p,k),
                 update.sigma = FALSE,update.sa = FALSE,update.w = FALSE,
                 maxiter = 20,tol = 0,drop.threshold = 0,verbose = FALSE)
b   <- with(out,rowSums(alpha * mu))

# Should be close to zero.
print(max(abs(fit$B - b)))

