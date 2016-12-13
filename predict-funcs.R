#' @title Predict method for vsm2s objects
#' 
#' @description Prediction with objects created by \code{\link{vsm2s}}
#' 
#' @param object object created by \code{\link{vsm2s}}
#' @param newx new predictor values (at which to predict)
#' @param which which of the fits to use for prediction (relevant if the \code{\link{vsm2s}} 
#' object includes multiple fits)
#' @param warn logical: Issue warning if \code{newx} truncated so as not to exceed range of original x?
#' @param ... other arguments
#' @importFrom fda eval.basis
#' 
#' @export
predict.vsm2s = function(object, newx, which=1, warn=TRUE, ...) {   
    if (warn & (min(newx)<min(object$x) | max(newx)>max(object$x)))
        warning("'newx' truncated so as not to exceed range of original x") 
    xtrunc = pmin(pmax(newx, min(object$x)), max(object$x))  
    eval.basis(xtrunc, object$basis.x) %*% object$xcoef[,,which]
}

#' @title Predict method for tps.fd objects
#' 
#' @description Prediction with objects created by \code{\link{tps.fd}}
#' 
#' @param object object created by \code{\link{tps.fd}}
#' @param newx new predictor values (at which to predict)
#' @param warn logical: Issue warning if \code{newx} truncated so as not to exceed range of original x?
#' @param ... other arguments
#' @importFrom fda eval.basis
#' 
#' @export
predict.tps.fd = function(object, newx, warn=TRUE, ...) {   
    if (warn & (min(newx)<min(object$x) | max(newx)>max(object$x)))
        warning("'newx' truncated so as not to exceed range of original x") 
    xtrunc = pmin(pmax(newx, min(object$x)), max(object$x))  
    tcrossprod(eval.basis(xtrunc, object$basis.x) %*% object$coef, object$B.f) 
}
