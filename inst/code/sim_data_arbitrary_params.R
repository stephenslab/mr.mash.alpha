###Set options
options(stringsAsFactors = F)

###Set seed
set.seed(123)

###Set parameters
n <- 600
p <- 100
p_causal <- 2
r <- 10

###Simulate V, B, X and Y
V <- mmbr:::create_cov_canonical(r, singletons=F, hetgrid=c(0.25))
V <- 0.8 * V[[1]]
Sigma <- mmbr:::create_cov_canonical(r, singletons=F, hetgrid=c(1))
scale_factor <- 0.8
B_causal <- MASS::mvrnorm(n=p_causal, mu=rep(0, r), Sigma=scale_factor*Sigma[[1]])
B <- matrix(0, ncol=r, nrow=p)
B[1:p_causal, ] <- B_causal
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
X <- scale(X, center=TRUE, scale=FALSE)
Y <- mr.mash.alpha:::sim_mvr(X, B, V)

###Build the mixture prior
grid <- seq(0.1, 1, 0.1)
S0 <- mr.mash.alpha:::compute_cov_canonical(ncol(Y), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid, zeromat=TRUE)
comps <- length(S0)
w0 <- rep(1/comps, comps)