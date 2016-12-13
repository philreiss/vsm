#' @title Tensor product smoothing for functional responses
#' 
#' @description Given a set of functional responses and a scalar predictor, this function implements
#' ordinary least squares and generalized least squares bivariate smooths, with respect to the 
#' function argument and the predictor
#'
#' @param Y \eqn{n \times L} matrix of functional responses
#' @param x \eqn{n}-dimensional predictor vector 
#' @param argvals function argument values 
#' @param k.f dimension of function-direction basis
#' @param k.x dimension of predictor-direction basis
#' @param norder.x order of splines in the predictor direction
#' @param basistype.f function-direction basis type for the functional direction
#' @param q.x,q.f penalty degree for the predictor- and function-direction bases
#' @param method "P-OLS" for penalized ordinary least squares, "P-GLS" for penalized 
#' generalized least squares
#' @param cov.method covariance estimation method for GLS: "mod.chol" for modified Choleski 
#' decomposition, or "naive" 
#' @param adapt logical: should the predictor-dimension smoothing parameter be chosen 
#' in a spatially adaptive manner? 
#' @param sp.f,sp.x smoothing parameter, passed to \code{\link[mgcv]{bam}}
#'
#' @return An object of class \code{tps.fd}, which is a list whose elements include
#' \item{coef}{spline coeffients}
#' \item{fitted}{fitted values}
#' \item{pwdf}{pointwise degrees of freedom}
#' 
#' @references Reiss, P. T., Huang, L., Chen, H., and Colcombe, S. Varying-smoother models for 
#' functional responses. Available at \url{https://arxiv.org/abs/1412.0778}
#'
#' Reiss, P. T., Huang, L., Wu, P.-S., Chen, H., and Colcombe, S. Pointwise influence matrices for
#' functional-response regression. Available at \url{https://works.bepress.com/phil_reiss/43/}
#'
#' @importFrom mgcv bam
#' @importFrom fda eval.basis eval.penalty
#' 
#' @export


tps.fd <- function(Y, x, argvals=NULL, k.f=25, 
                 k.x=15, norder.x=4, basistype.f = "bspline", q.f=2, 
                 q.x=norder.x-2, method="P-OLS", cov.method="mod.chol", adapt=FALSE, 
                 sp.f=-1, sp.x=-1) {
    if (adapt & basistype.f!="bspline") stop("Adaptive smoothing is implemented only for basistype.f=='bspline'.")
    if (adapt && k.f/2==floor(k.f/2)) warning("'k.f' should ideally be odd for adaptive smoothing.")
    if (adapt && q.f!=2) warning("'q.f' should ideally be 2 for adaptive smoothing.")
    if ((sp.f>0|sp.x>0) & (method!="P-OLS" | adapt)) warning("Specified smoothing parameters used only for preliminary OLS fit") #3/5
  
    if (is.null(argvals)) argvals = 1:ncol(Y)
    n = nrow(Y)
    L = ncol(Y)
    Init = vsmInit(Y=Y, x=x, argvals=argvals, basistype.f=basistype.f, 
                   k.f=k.f, q.f=q.f, norder.x=norder.x, k.x=k.x, q.x=q.x)
    basis.f=Init$basis.f; basis.x=Init$basis.x;
    B.f=Init$B.f; P.f=Init$P.f; B.x=Init$B.x; P.x=Init$P.x; 
    B.der.x <- eval.basis(x, basis.x, Lfdobj=1)
  
    # For tensor product penalty
    M.x <- eval.penalty(basis.x, Lfdobj=0)
    M.f <- eval.penalty(basis.f, Lfdobj=0)
    Pf.Mx <- P.f %x% M.x
    Mf.Px <- M.f %x% P.x
    x.design <- B.f %x% B.x
    penlist = list(Pf.Mx)
    if (!adapt) penlist[[2]] = Mf.Px
    else {
  	    lilbasis = create.bspline.basis(rangeval=basis.f$rangeval, breaks=sort(c(basis.f$rangeval, median(basis.f$params))))
  	    penmats = multipen.nc(basis.f, lilbasis)
        for (l in 1:5) penlist[[l+1]] = penmats[[l]] %x% P.x
    }
    if (!adapt)  mod.ols <- bam(as.vector(Y)~x.design-1, paraPen=list(x.design=penlist, sp=c(sp.f, sp.x)),method="REML") #3/5
    else mod.ols <- bam(as.vector(Y)~x.design-1, paraPen=list(x.design=penlist),method="REML") #3/5
    sp = mod.ols$sp
  
    if (method=="P-OLS"){
        tpcoef <- matrix(mod.ols$coef, k.x, k.f)
        fitder = tcrossprod(B.der.x %*% tpcoef, B.f) 
        xcoef <- tcrossprod(tpcoef, B.f)   # get roughness
        roughness = colSums(xcoef * (P.x %*% xcoef))
        penpart = 0
        for (r in 1:length(penlist)) penpart = penpart + mod.ols$sp[r] * penlist[[r]]
        W <- (matrix(colSums(B.f),1)%x%B.x) %*% solve(crossprod(B.f)%x%crossprod(B.x)+penpart)
        pwdf <- rowSums((diag(L)%x%matrix(1,1,n)) %*% ((B.f%x%B.x)*(rep(1,L)%x%W)))
    }
  
    fitted <- matrix(mod.ols$fit, n, L)
  
    if (method=="P-GLS"){
        residmat <- Y-fitted
        ### estimating the covariance matrix for prewhitening
        if (cov.method == "mod.chol") {
            p = ncol(residmat)
            res.cent <- scale(residmat, TRUE, FALSE)
            sqrt.prec.list = list()
            lwstat = lwpval = c()
            for (nband in 1:(min(n-2,p-1))) {
                TT = diag(p)
                Ddiag = rep(0, p)
                Ddiag[1] = var(res.cent[, 1])
                for (k in 2:p) {
                    qrResCent <- qr(res.cent[, max(1, k - nband):(k - 1)])
                    TT[k, max(1, k - nband):(k - 1)] <- (-qr.coef(qrResCent, res.cent[, k]))
                    Ddiag[k] <- var(qr.resid(qrResCent, res.cent[, k]))
                }
                prec <- scale(t(TT), FALSE, Ddiag) %*% TT
                svdp <- eigen(prec, symmetric = TRUE)
                sqrt.prec.list[[nband]] = svdp$vectors %*% 
                                             tcrossprod(diag(sqrt(svdp$values)), svdp$vectors)
                lwprec = lw.test(residmat %*% sqrt.prec.list[[nband]])
                lwstat[nband] = lwprec$stat
                lwpval[nband] = lwprec$pvalue
                if (lwstat[nband] < -5) break
            }
            nband.best = which.max(lwpval)
            cat("Using half-bandwidth", nband.best, "for precision matrix of residuals\n")
            sqrt.prec <- sqrt.prec.list[[nband.best]]
        }
        else if (cov.method == "naive") {
            if (nrow(residmat) < ncol(residmat)) 
            stop("Sample covariance matrix of residuals is singular.")
            svd.cov.mle <- svd(cov(residmat) * (n - 1)/n)
            sqrt.prec <- tcrossprod(scale(svd.cov.mle$u, FALSE, sqrt(svd.cov.mle$d)), svd.cov.mle$u)
        }

        x.design.prewhiten <- (sqrt.prec%*%B.f) %x% B.x
        mod.gls <- bam(as.vector(Y %*% sqrt.prec)~x.design.prewhiten-1,  
                       paraPen=list(x.design.prewhiten=penlist), method="REML")
        tpcoef <- matrix(mod.gls$coef, k.x, k.f)
        fitted <- tcrossprod(B.x %*% tpcoef, B.f)  
        fitder = tcrossprod(B.der.x %*% tpcoef, B.f)
        sp = mod.gls$sp
    
        xcoef <- tcrossprod(tpcoef, B.f)   
        roughness = colSums(xcoef * (P.x %*% xcoef))
    
        penpart = 0
        for (r in 1:length(penlist)) penpart = penpart + mod.gls$sp[r] * penlist[[r]]
        W <- (matrix(colSums(crossprod(sqrt.prec)%*%B.f),1)%x%B.x) %*%
                     solve(crossprod(sqrt.prec%*%B.f)%x%crossprod(B.x)+penpart)
        pwdf <- rowSums((diag(L)%x%matrix(1,1,n)) %*% ((B.f%x%B.x)*(rep(1,L)%x%W)))
    }
  
    fitt <- list(coef=tpcoef, fitted=fitted, pwdf = pwdf, x = x, basis.x = basis.x, 
               basis.f = basis.f, B.x = B.x, B.f = B.f, P.x = P.x, B.der.x = B.der.x,
               P.f = P.f, argvals = argvals, Y = Y, fitder=fitder, roughness=roughness, sp=sp)
    class(fitt) <- "tps.fd"
    fitt
}
