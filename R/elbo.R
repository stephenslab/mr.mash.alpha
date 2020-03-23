###Compute ELBO from intermediate components
compute_ELBO_fun <- function(Rbar, V, Vinv, ldetV, var_part_tr_wERSS, neg_KL){
  n <- nrow(Rbar)
  r <- ncol(Rbar)
  # tr_wERSS <- tr(Vinv%*%(crossprod(Rbar))) + var_part_tr_wERSS
  tr_wERSS <- sum(Vinv*(crossprod(Rbar))) + var_part_tr_wERSS
  ELBO <- -log(n)/2 - (n*r)/2*log(2*pi) - n/2 * ldetV - 0.5*tr_wERSS + neg_KL
  
  return(ELBO)
}

###Compute intermediate components of the ELBO
compute_ELBO_terms <- function(var_part_tr_wERSS, neg_KL, x_j, Rbar_j, bfit, xtx, Vinv){
  mu1_mat <- matrix(bfit$mu1, ncol=1)
  # var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%bfit$S1)*xtx)
  # neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, Rbar_j)%*%tcrossprod(matrix(x_j, ncol=1), mu1_mat))+
  #                                        tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*xtx))
  ##Equivalent to the above but more efficient
  var_part_tr_wERSS <- var_part_tr_wERSS + (sum(Vinv*bfit$S1)*xtx)
  neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*sum(tcrossprod(Vinv, Rbar_j)*t(tcrossprod(matrix(x_j, ncol=1), mu1_mat)))+
                                         sum(Vinv*(bfit$S1+tcrossprod(mu1_mat)))*xtx))
  
  
  return(list(var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
}

