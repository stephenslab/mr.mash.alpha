#' @title Predict future observations from mr.mash fit.
#' 
#' @param object a mr.mash fit.
#' 
#' @param newx a new value for X at which to do predictions.
#'
#' @param \dots Additional arguments (not used).
#' 
#' @return Matrix of predicted values.
#'
#' @importFrom stats predict
#'
#' @method predict mr.mash
#' @export
#' @export predict.mr.mash
#' 
predict.mr.mash <- function(object, newx, ...){
  if(!is.matrix(newx))
    stop("X must be a matrix.")
  if (any(is.na(newx)))
    stop("X must not contain missing values.")
  return(with(object,addtocols(newx %*% mu1,intercept)))
}

#' @title Extract coefficients from mr.mash fit.
#' 
#' @param object a mr.mash fit.
#'
#' @param \dots Other arguments (not used).
#' 
#' @return (p+1) x r matrix of coefficients.
#' 
#' @importFrom stats coef
#' 
#' @method coef mr.mash
#' @export coef.mr.mash
#' @export
#' 
coef.mr.mash <- function(object, ...){
  coeffs <- rbind(object$intercept, object$mu1)
  rownames(coeffs)[1] <- "(Intercept)"
  return(coeffs)
}
