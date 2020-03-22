###Compute ELBO from intermediate components
compute_ELBO_fun <- function(rbar, V, Vinv, ldetV, var_part_tr_wERSS, neg_KL){
  n <- nrow(rbar)
  R <- ncol(rbar)
  # tr_wERSS <- tr(Vinv%*%(crossprod(rbar))) + var_part_tr_wERSS
  tr_wERSS <- sum(Vinv*(crossprod(rbar))) + var_part_tr_wERSS
  ELBO <- -log(n)/2 - (n*R)/2*log(2*pi) - n/2 * ldetV - 0.5*tr_wERSS + neg_KL
  
  return(ELBO)
}

###Compute intermediate components of the ELBO
compute_ELBO_terms <- function(var_part_tr_wERSS, neg_KL, x_j, rbar_j, bfit, xtx, Vinv){
  mu1_mat <- matrix(bfit$mu1, ncol=1)
  # var_part_tr_wERSS <- var_part_tr_wERSS + (tr(Vinv%*%bfit$S1)*xtx)
  # neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*tr(tcrossprod(Vinv, rbar_j)%*%tcrossprod(matrix(x_j, ncol=1), mu1_mat))+
  #                                        tr(Vinv%*%(bfit$S1+tcrossprod(mu1_mat)))*xtx))
  ##Equivalent to the above but more efficient
  var_part_tr_wERSS <- var_part_tr_wERSS + (sum(Vinv*bfit$S1)*xtx)
  neg_KL <- neg_KL + (bfit$logbf +0.5*(-2*sum(tcrossprod(Vinv, rbar_j)*t(tcrossprod(matrix(x_j, ncol=1), mu1_mat)))+
                                         sum(Vinv*(bfit$S1+tcrossprod(mu1_mat)))*xtx))
  
  
  return(list(var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
}

###Compute variance part of the ERSS
compute_var_part_ERSS <- function(var_part_ERSS, bfit, xtx){
  var_part_ERSS <- var_part_ERSS + (bfit$S1*xtx)
  
  return(var_part_ERSS)
}

