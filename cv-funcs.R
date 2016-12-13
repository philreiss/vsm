#' @title Cross-validation for two-step smoothing
#' 
#' @description Cross-validation function for \code{\link{vsm2s}}. Useful for choosing optimal 
#' function-direction smoothing parameter.
#'
#' @param Y \eqn{n \times L} matrix of functional responses
#' @param x \eqn{n}-dimensional predictor vector 
#' @param lsp.f candidate log smoothing parameter in the function direction
#' @param nfpc number of functional principal components
#' @param n.folds number of folds
#' @param seed seed for random number generation
#' @param argvals function argument values
#' @param raw use raw (first-step) coefficients
#' @param ... other arguments passed to \code{\link{vsm2s}}
#' @importFrom fda eval.basis
#' 
#' @export
cv.vsm2s = function(Y, x, lsp.f=NULL, nfpc=NULL, n.folds = 5, seed = NULL, argvals=NULL, raw = FALSE, ...) { 
	if (!is.null(seed)) set.seed(seed)   
    min.folds = sample(order(x)[1:n.folds])
    max.folds = sample(order(x)[(length(x)-n.folds+1):length(x)])
    rest.folds = split(sample((1:length(x))[-c(min.folds, max.folds)]), rep(1:n.folds, length=length(x)-2*n.folds))
    if (!is.null(lsp.f)) nfit = length(lsp.f)
    else if (!is.null(nfpc)) nfit = length(nfpc)
    else nfit = 1
    cvarray = array(dim = c(ncol(Y), nfit, n.folds))
    if (!is.null(lsp.f)) dimnames(cvarray) = list(NULL, lsp.f, 1:n.folds)
    else if (!is.null(nfpc)) dimnames(cvarray) = list(NULL, nfpc, 1:n.folds)
    
    for (k in 1:n.folds)    {
        cat("####### Fold",k,"#######\n")
        val = c(min.folds[k], rest.folds[[k]], max.folds[k])
        trainmod = vsm2s(Y=Y[-val, ], x=x[-val], argvals=argvals, lsp.f = lsp.f, nfpc=nfpc, ...) 
        if (raw) {
            for (j in 1:nfit) {
                xtrunc = pmin(pmax(x[val], min(trainmod$x)), max(trainmod$x))  
                temp = eval.basis(xtrunc, trainmod$basis.x)
                cvarray[,j,k] = colSums((Y[val, ] - temp%*%trainmod$coef.raw)^2)
            }
        }
        else {
            for (j in 1:nfit)
                cvarray[,j,k] = colSums((Y[val, ] - predict(trainmod, as.vector(x[val]), which=j, warn=FALSE))^2)
        }
    }
    
    cvmat = sqrt(apply(cvarray,1:2,sum) / length(x))
    cvvec = sqrt(rowSums(apply(cvarray, 2:3, sum)) / (length(Y)))
    
    if (length(lsp.f)>1) {
        if (which.min(cvvec)==which.min(lsp.f)) warning("CV minimized at lowest smoothing parameter value")
        if (which.min(cvvec)==which.max(lsp.f)) warning("CV minimized at highest smoothing parameter value")
    }
    if (length(nfpc)>1) {
        if (which.min(cvvec)==which.min(nfpc)) warning("CV minimized at lowest number of fpc")
        if (which.min(cvvec)==which.max(nfpc)) warning("CV minimized at highest number of fpc")
    }
    
    # obtain the optimal idx of the candidate smoothing parameter
    which.sp = which.min(apply(cvmat, 2,sum))
    # obtain the optimal idx vec pointwisely
    pw.which.sp = apply(cvmat, 1, which.min)
  
    return(list(cvarray = cvarray, cvmat = cvmat, cvvec = cvvec, which.sp = which.sp, pw.which.sp = pw.which.sp))
}
 
#' @title Cross-validation for tensor product smoothing for functional responses
#' 
#' @description Cross-validation function for \code{\link{tps.fd}}. Can be used for choosing optimal 
#' function-direction smoothing parameter.
#'
#' @param Y \eqn{n \times L} matrix of functional responses
#' @param x \eqn{n}-dimensional predictor vector 
#' @param n.folds number of folds
#' @param seed seed for random number generation
#' @param argvals function argument values
#' @param ... other arguments passed to \code{\link{tps.fd}}
#' 
#' @export
cv.tps.fd = function(Y, x, n.folds = 5, seed = NULL, argvals=NULL, ...) { 
	# Can this whole function be merged with cv.vsm2s?
	if (!is.null(seed)) set.seed(seed)   
    min.folds = sample(order(x)[1:n.folds])
    max.folds = sample(order(x)[(length(x)-n.folds+1):length(x)])
    rest.folds = split(sample((1:length(x))[-c(min.folds, max.folds)]), rep(1:n.folds, length=length(x)-2*n.folds))
    nfit = 1
    cvarray = array(dim = c(ncol(Y), nfit, n.folds))
    
    for (k in 1:n.folds)    {
        cat("####### Fold",k,"#######\n")
        val = c(min.folds[k], rest.folds[[k]], max.folds[k])
        trainmod = tps.fd(Y=Y[-val, ], x=x[-val], argvals=argvals, ...) 
        for (j in 1:nfit)
            cvarray[,j,k] = colSums((Y[val, ] - predict(trainmod, as.vector(x[val]), warn=FALSE))^2)            
    }
    
    cvmat = sqrt(apply(cvarray,1:2,sum) / length(x))
    cvvec = sqrt(rowSums(apply(cvarray, 2:3, sum)) / (length(Y)))
  
    return(list(cvarray = cvarray, cvmat = cvmat, cvvec = cvvec))
}
