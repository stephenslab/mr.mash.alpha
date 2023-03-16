#' @title Predict future observations from mr.mash.rss fit.
#' 
#' @param object a mr.mash.rss fit.
#' 
#' @param newx a new value for X at which to do predictions.
#'
#' @param \dots Additional arguments (not used).
#' 
#' @return Matrix of predicted values.
#'
#' @export
#' @export predict.mr.mash.rss
#' 
predict.mr.mash.rss <- function(object, newx, ...){
  if(!is.matrix(newx))
    stop("X must be a matrix.")
  if (any(is.na(newx)))
    stop("X must not contain missing values.")
  
  if(!is.na(object$intercept))
    return(with(object,addtocols(newx %*% mu1,intercept)))
  else
    return(with(object,newx %*% mu1))
}

#' @title Extract coefficients from mr.mash.rss fit.
#' 
#' @param object a mr.mash fit.
#'
#' @param \dots Other arguments (not used).
#' 
#' @return p x r or (p+1) x r matrix of coefficients,
#'  depending on whether an intercept was computed.
#' 
#' @export
#' @export coef.mr.mash.rss
#' 
coef.mr.mash.rss <- function(object, ...){
  if(!is.na(object$intercept)){
    coeffs <- rbind(object$intercept, object$mu1)
    rownames(coeffs)[1] <- "(Intercept)"
  } else {
    coeffs <- object$mu1
  }
  
  return(coeffs)
}
