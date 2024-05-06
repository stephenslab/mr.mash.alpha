###Compute ELBO from intermediate components
compute_ELBO_rss_fun <- function(n, RbartRbar, Vinv, ldetV, var_part_tr_wERSS, neg_KL){
  r <- ncol(RbartRbar)
  
  tr_wERSS <- sum(Vinv*RbartRbar) + var_part_tr_wERSS

  ELBO <- -log(n)/2 - (n*r)/2*log(2*pi) - n/2 * ldetV - 0.5*(tr_wERSS) + neg_KL
  
  return(ELBO)
}

###Compute intermediate components of the ELBO
compute_ELBO_rss_terms <- function(var_part_tr_wERSS, neg_KL, xtRbar_j, bfit, xtx, Vinv){
  mu1_mat <- matrix(bfit$mu1, ncol=1)

  var_part_tr_wERSS <- var_part_tr_wERSS + (sum(Vinv*bfit$S1)*xtx)
  
  mu1xtRbar_j <- mu1_mat%*%xtRbar_j
  
  Cm <- - mu1xtRbar_j - t(mu1xtRbar_j) + tcrossprod(mu1_mat)*xtx + bfit$S1*xtx
  
  neg_KL <- neg_KL + (bfit$logbf + 0.5*(sum(Vinv*Cm)))
  
  return(list(var_part_tr_wERSS=var_part_tr_wERSS, neg_KL=neg_KL))
}

