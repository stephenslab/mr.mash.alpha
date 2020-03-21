#' @title Predict future observations from mr.mash fit.
#' 
#' @param object a mr.mash fit.
#' 
#' @param newx a new value for X at which to do predictions.
#' 
#' @return a matrix of predicted values.
#' 
#' @export
#' 
predict.mr.mash <- function(object, newx){
  if(!is.matrix(newx))
    stop("X must be a matrix.")
  if (any(is.na(newx)))
    stop("X must not contain missing values.")
  return(with(object,addtocols(newx %*% mu1,intercept)))
}
