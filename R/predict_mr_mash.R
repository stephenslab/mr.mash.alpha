#' @title Predict future observations from mr.mash fit.
#' @param object a mr.mash fit.
#' @param newx a new value for X at which to do predictions.
#' @return a matrix of predicted values.
#' @S3method predict mr.mash
#' @export
#' 
predict.mr.mash <- function(object, newx){
  if(!is.matrix(newx)){
    stop("X must be a matrix.")
  }
  if (any(is.na(newx))) {
    stop("X must not contain missing values.")
  }
  if(!is(object,"mr.mash")){
    stop("Input argument object must be an instance of class \"mr.mash\".")
  }
  
  Yhat <- matrix(rep(object$intercept, each=nrow(newx)), ncol=ncol(object$mu1)) + newx%*%object$mu1
  
  return(Yhat)
}
