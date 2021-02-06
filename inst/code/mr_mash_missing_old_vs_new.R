# An illustration of the mr_mash_simple implementation applied to a
# small, simulated data set.
suppressMessages(library(MBSP))
source("../../R/misc.R")
source("../../R/mr_mash_simple.R")

# SCRIPT PARAMETERS
# -----------------
# Number of samples (n) and number of predictors (p).
n <- 50
p <- 100

# Residual covariance matrix.
V <- rbind(c(1.0,0.2),
           c(0.2,0.4))
r <- nrow(V)

# True effects used to simulate the data.
B <- rbind(c(-2.0, -1.5),
           c( 1.0,  1.0),
           matrix(0,p - 2,r))

# Covariances in the mixture-of-normals prior on the regression
# coefficients.
S0 <- list(k1 = rbind(c(3,0),
                      c(0,3)),
           k2 = rbind(c(4,2),
                      c(2,4)),
           k3 = rbind(c(6,3.5),
                      c(3.5,4)),
           k4 = rbind(c(5,0),
                      c(0,0)),
           k5 = rbind(c(0,0),
                      c(0,0)))

# The mixture weights in the mixture-of-normals prior on the
# regression coefficients.
w0 <- c(0.1,0.2,0.2,0.1,0.4)

# SIMULATE DATA
# -------------
set.seed(1)
X <- matrix(rnorm(n*p),n,p)
X <- scale(X,scale = FALSE)

# Simulate Y ~ MN(X*B,I,V). Note that matrix.normal from the MBSP
# package appears to be much faster than rmatrixnorm from the
# MixMatrix package.
Y <- matrix.normal(X %*% B,diag(n),V)

# Assign some missing values in Y
Y_miss <- Y
Y_miss[1:5, 1] <- NA
Y_miss[11:15, 2] <- NA

# FIT MR-MASH MODEL ALLOWING FOR MISSING Ys
# -----------------
# Run 20 co-ordinate ascent updates (older function).
B0  <- matrix(0,p,r)
fit_miss <- mr_mash_simple_missing_Y(X,Y_miss,V,S0,w0,B0,20,update_w0 = TRUE,
                                     update_V = TRUE,verbose = TRUE)

# Run 20 co-ordinate ascent updates (newer function).
fit_miss_1 <- mr_mash_simple_missing_Y_1(X,Y_miss,V,S0,w0,B0,20,update_w0 = TRUE,
                                         update_V = TRUE,verbose = TRUE)


# Compare the posterior mean estimates of the regression coefficients
plot(fit_miss$B,fit_miss_1$B,pch = 20,xlab = "Old",ylab = "New", main = "Effects")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

# Compare the estimate of the residual covariance
plot(fit_miss$V,fit_miss_1$V,pch = 20,xlab = "Old",ylab = "New", main = "Residual covariance")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

# Compare the estimate of the prior weights
plot(fit_miss$w0,fit_miss_1$w0,pch = 20,xlab = "Old",ylab = "New", main = "Prior weights")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")

# Compare the estimate of the ELBO
plot(fit_miss$ELBO,fit_miss_1$ELBO,pch = 20,xlab = "Old",ylab = "New", main = "ELBO")
abline(a = 0,b = 1,col = "skyblue",lty = "dotted")



