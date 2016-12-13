#' @title Varying-smoother model via two-step smoothing
#' 
#' @description Fits a nonlinear function-on-scalar regression model by (1) smoothing in the 
#' "predictor direction" at each point along the function, then (2) smoothing the result in the 
#' "function direction" (along the function). Supports only a single predictor.
#'
#' @param Y \eqn{n \times L} matrix of functional responses
#' @param x \eqn{n}-dimensional predictor vector 
#' @param argvals function argument values 
#' @param k.f dimension of function-direction basis
#' @param k.x dimension of predictor-direction basis
#' @param lsp.x,lsp.f candidate log smoothing parameters in the predictor and function direction
#' @param norder.x order of splines in the predictor direction
#' @param basistype.f function-direction basis type for the functional direction
#' @param q.x,q.f penalty degree for the predictor- and function-direction bases
#' @param varest logical: should variance be estimated?
#' @param nfpc number of functional principal components
#' @param nfpc.default default number of functional principal components
#' @param pve.cutoff proportion of variance explained; used to choose number of FPCs
#'
#' @return An object of class \code{vsm2s}, which is a list whose elements include
#' \item{tpcoef,xcoef,coef.raw,coef.smoothed}{several types of spline coeffients}
#' \item{fitted}{fitted values}
#' \item{pwdf}{pointwise degrees of freedom}
#'
#' @references Reiss, P. T., Huang, L., Chen, H., and Colcombe, S. Varying-smoother models for 
#' functional responses. Available at \url{https://arxiv.org/abs/1412.0778}
#'
#' Reiss, P. T., Huang, L., Wu, P.-S., Chen, H., and Colcombe, S. Pointwise influence matrices for
#' functional-response regression. Available at \url{https://works.bepress.com/phil_reiss/43/}
#'
#' @importFrom fda pca.fd Data2fd
#' 
#' @export

vsm2s = function(Y, x, argvals=NULL, k.f=30+(basistype.f=="fourier"), k.x=15, 
                 lsp.x=seq(-10, 20, length.out=100), lsp.f=NULL, norder.x=4, basistype.f = "bspline", 
                 q.f=2, q.x=norder.x-2, varest=FALSE, nfpc=NULL, nfpc.default=10, 
                 pve.cutoff=.99) {

    if (is.null(argvals)) argvals = 1:ncol(Y)
    n = nrow(Y)
    L = ncol(Y)
    Init = vsmInit(Y=Y, x=x, argvals=argvals, basistype.f=basistype.f, k.f=k.f, q.f=q.f, norder.x=norder.x, k.x=k.x, q.x=q.x)
    basis.f=Init$basis.f; basis.x=Init$basis.x
    B.f=Init$B.f; P.f=Init$P.f; B.x=Init$B.x; P.x=Init$P.x 
    RinvU.x=Init$RinvU.x; A.x=Init$A.x; tau=Init$tau
    B.der.x <- eval.basis(x, basis.x, Lfdobj=1)
    
    if (is.null(lsp.f) && is.null(nfpc)) method = "smooth"
    else if (is.null(lsp.f) && !is.null(nfpc)) method = "FPC"
    else if (!is.null(lsp.f) && is.null(nfpc)) method = "smooth"
    else if (!is.null(lsp.f) && length(nfpc)==1) method = "penFPC"
    else stop("The combination of lsp.f and nfpc is not allowed")
   
    # Step 1: Massively parallel smoothing        
    svdP = svd(P.x)
    Us = svdP$u[ , 1:(k.x-q.x)]
    Un = if (q.x>0) svdP$u[ , -(1:(k.x-q.x))] else rep(0,k.x)
    X = B.x %*% Un
    Z = scale(B.x %*% Us, FALSE, sqrt(svdP$d[1:(k.x-q.x)]))  # see Wood (2004)
    svdZZ = svd(tcrossprod(Z))
  	d = svdZZ$d; d[-(1:ncol(Z))] = 0
  	X. = crossprod(svdZZ$u, X)
  	X.X. = crossprod(X.)
  	Y. = crossprod(svdZZ$u, Y)
  	m1 = Y. - X. %*% solve(X.X., crossprod(X., Y.))
  	rll = Vectorize(function(log.sp) {
    		  ev.v = 1 + exp(-log.sp)*d
    		  di = diag(1/sqrt(ev.v))
    		  X.. = di %*% X.
    		  X..X.. = crossprod(X..)
    		  Y.. = di %*% Y.
    		  m2 = Y.. - X.. %*% solve(X..X.., crossprod(X.., Y..))
    		    - (n-ncol(X)) * log(colSums(Y.. * m2)) - sum(log(ev.v)) - log(det(X..X..))
	      })
	tabl = rll(lsp.x)
  	colnames(tabl) = lsp.x
  	best.lsp.x = lsp.x[as.numeric(apply(tabl, 1, which.max))]
    M = 1 / (1 + tau %o% exp(best.lsp.x))
    coef.raw = RinvU.x %*% (M * crossprod(A.x, Y))
    coef.smoothed = M.smoothed = NULL
    pwdf.raw = colSums(M)
    
    ########################################################

    # Step 2
    nfit = if (is.null(lsp.f)) max(length(nfpc),1) else max(length(lsp.f),1)
    tpcoef = array(dim=c(k.x, k.f, nfit))
    xcoef = array(dim=c(k.x, L, nfit))
    fitted = fitder = array(dim=c(n, L, nfit))
    pwdf = matrix(NA, L, nfit) 

    if (method=="smooth") {   
        varprop = NULL
        
        # Demmler-Reinsch in function direction
        Rinv.f = solve(chol(crossprod(B.f)))
        svdRPR.f = svd(crossprod(Rinv.f, P.f %*% Rinv.f))
        RinvU.f = Rinv.f %*% svdRPR.f$u; ss.f = svdRPR.f$d
        A.f = B.f %*% RinvU.f
        leftmat = coef.raw %*% A.f
        for (k in 1:nfit) {
   	        skale =  1 + exp(lsp.f[k]) * ss.f
       	    tpcoef[ , , k] = tcrossprod(scale(leftmat, FALSE, skale), RinvU.f)
            S.f = tcrossprod(scale(A.f, FALSE, skale), A.f)
   	        pwdf[ , k] = crossprod(pwdf.raw, S.f)
       	}
    } 
      
    else if (method %in% c("FPC", "penFPC")) {   
    		meanmat = outer(rep(1, k.x), colMeans(Y))
        tpcoef.meanpart = meanmat %*% B.f %*% solve(crossprod(B.f))
    		if (method=="FPC") { 
    	    		fpcobj = pca.fd(Data2fd(argvals=argvals, y=t(Y), basisobj=basis.f), nharm=max(nfpc))
        		varprop = fpcobj$varprop
    	    		V.all = fpcobj$harmonics$coef
            BV.all = B.f %*% V.all
            for (k in 1:nfit) {
            		V = V.all[ , 1:nfpc[k]]
    		    		BV = BV.all[ , 1:nfpc[k]]
    		    		tpcoef[ , , k] = tpcoef.meanpart + (coef.raw - meanmat) %*% BV %*% solve(crossprod(BV), t(V)) 
    		    		pwdf[ , k] = 1 + BV %*% solve(crossprod(BV), t(BV)) %*% (pwdf.raw - 1)
            }
    		}
        else if (method=="penFPC") { 
    	    		fpcobj = pca.fd(Data2fd(argvals=argvals, y=t(Y), basisobj=basis.f), nharm=if (nfpc>0) nfpc else min(n-1, k.f))
    	    		varprop = fpcobj$varprop
    	    		nfpc.use = if (nfpc>0) nfpc else max(nfpc.default, which.min(cumsum(varprop) > pve.cutoff))
    	    		cat("penFPC method uses", nfpc.use, "components\n")
    	    		V = fpcobj$harmonics$coef[ , 1:nfpc.use]
            BV = B.f %*% V
            for (k in 1:nfit) {
    	            tpcoef[ , , k] = tpcoef.meanpart + (coef.raw - meanmat) %*% BV %*% solve(crossprod(BV) + exp(lsp.f[k])*crossprod(V, P.f %*% V), t(V)) 
    		        pwdf[ , k] = 1 + BV %*% solve(crossprod(BV)+exp(lsp.f[k])*crossprod(V, P.f %*% V), t(BV)) %*% (pwdf.raw - 1)
            }
        }
    }   
    
    roughness <- matrix(NA, L, nfit)   
    
    for (k in 1:nfit) {
        xcoef[ , , k] <- if (is.null(coef.smoothed)) tcrossprod(tpcoef[ , ,k], B.f) else coef.smoothed
        fitted[ , , k] = B.x %*% xcoef[ , , k]
        fitder[ , , k] = B.der.x %*% xcoef[ , , k]
        roughness[ ,k] = colSums(xcoef[ , , k] * (P.x %*% xcoef[ , , k]))      
    }
    
    pwgcv = apply(fitted, 3, function(x) colMeans((Y - x)^2))   / (1 - pwdf/n)^2
  
   	fitt = list(tpcoef = tpcoef, xcoef=xcoef, fitted = fitted, fitder = fitder,
    	        		pwdf.raw = pwdf.raw, pwdf = pwdf, pwlsp = best.lsp.x,
   	            coef.raw=coef.raw, coef.smoothed=coef.smoothed,
   	            x = x, basis.x = basis.x, M=M, M.smoothed=M.smoothed,
   	            basis.f = basis.f, B.x = B.x, B.f = B.f, P.x = P.x, B.der.x = B.der.x,
   	            P.f = P.f, argvals = argvals, Y = Y, varprop = varprop, 
                roughness = roughness, pwgcv = pwgcv)

    if (varest) {
    		if (nfit>1) stop("Variance estimation not currently implemented with multiple candidate fits.")
    		if (method!="smooth") stop("Variance estimation currently implemented only for B-spline smoothing method.")

    		# Estimate residual var-cov matrix
    		cat("Estimating residual variance...\n")
    		res = Y - fitted[,,1]
    		sigmahat = cov(res)

    		cat("Extracting matrix for pointwise CIs...\n")
    		svdsig = svd(sigmahat)
    		halfsig = diag(sqrt(pmax(svdsig$d, 0))) %*% svdsig$v
    		# Can the next line be done more efficiently using A.f, etc.?
    	# T2 = if (is.null(lsp.f)) SginvB.f else B.f %*% solve(crossprod(B.f)+exp(lsp.f)*P.f)
    	# See old versions of code (e.g., frs.R) for definition of SginvB.f
    	if (is.null(lsp.f)) stop("Variance estimation not implemented if 'lsp.f' is NULL")
    	T2 = B.f %*% solve(crossprod(B.f)+exp(lsp.f)*P.f)
		fitt$TT = scale(halfsig %x% diag(k.x), FALSE, 1/as.vector(M)) %*% (T2 %x% t(RinvU.x))
   	}
   	    
    class(fitt) = "vsm2s"
    fitt
}
