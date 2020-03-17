# Small script to illustrate how computation of posterior quantities
# and Bayes factors can be made more efficient by pre-computing some
# things in advance.
library(mvtnorm)

# Returns the log-density of the multivariate normal with zero mean
# and covarirance S at x.
ldmvnorm <- function (x, S)
  dmvnorm(x,sigma = S,log = TRUE)

# Simulate data.
set.seed(1)
n <- 20
x <- rnorm(n)
Y <- matrix(rnorm(2*n),n,2)
V <- rbind(c(1.0,0.2),
           c(0.2,0.4))
S0 <- rbind(c(6,3.5),
            c(3.5,4))

# Least-squares calculations.
xx   <- sum(x^2)
bhat <- drop(x %*% Y)/xx
S    <- V/xx

# Bayesian calculations.
I     <- diag(2)
S1    <- S0 %*% solve(I + solve(S) %*% S0)
mu1   <- drop(S1 %*% solve(S,bhat))
logbf <- ldmvnorm(bhat,S0 + S) - ldmvnorm(bhat,S)

# Store the Bayesian calculations in a list.
dat <- list(mu1   = mu1,
            S1    = S1,
            logbf = logbf)

# Pre-calculations (that don't depend on X). Here, U0 is the prior
# covariance of the "transformed" data.
P   <- solve(V)
R   <- chol(V)
U0  <- solve(t(R)) %*% S0 %*% solve(R)
out <- eigen(U0)
d   <- out$values
Q   <- out$vectors
QR  <- t(Q) %*% R

# Same Bayesian calculations as before, but done more "efficiently" by
# making use of the "pre-calculations". Here, U1 is the posterior
# covariance of the "transformed" data.
D   <- diag(d/(1 + xx*d))
dx  <- d/(1 + xx*d)
B   <- sqrt(dx) * QR
S1  <- crossprod(B) # <-- Most expensive operation.
mu1 <- drop(B %*% (t(B) %*% (P %*% (bhat/xx))))

# Compare the two calculations.
print(mu1 - dat$mu1)
print(S1 - dat$S1)
