#' @importFrom fda create.fourier.basis create.bspline.basis eval.basis getbasispenalty
#' @importFrom mgcv gam
vsmInit = function(Y, x, argvals, basistype.f, k.f, q.f, norder.x, k.x, q.x)  {

    # Function-space basis and penalty
    if (basistype.f %in% c('fourier', 'bspline')) {
    	if (basistype.f=='fourier') basis.f = create.fourier.basis(range(argvals), k.f)
    	else if (basistype.f=='bspline') basis.f = create.bspline.basis(range(argvals), k.f)
    	B.f = eval.basis(argvals, basis.f)
    	P.f = getbasispenalty(basis.f, q.f)
    }
    else {     # use gam(..., fit=FALSE)
    	       # for some reason, basistype.f="cc" doesn't work
        unfit = gam(Y[1, ] ~ s(argvals, bs=basistype.f, k=k.f), fit=FALSE)
        B.f = unfit$X
        P.f = matrix(0, k.f, k.f)
        begin = unfit$off
        end = begin + nrow(unfit$S[[1]]) - 1
        P.f[begin:end, begin:end] = unfit$S[[1]]
        q.f = unfit$smooth[[1]]$null.space.dim
        basis.f = NULL
    }

    # Predictor-space basis and penalty
    basis.x = create.bspline.basis(range(x), breaks = quantile(x, seq(0,1,length.out=k.x-norder.x+2)), norder=norder.x)
    B.x = eval.basis(x, basis.x)
    P.x = getbasispenalty(basis.x, q.x)
	  
    # Demmler-Reinsch in predictor (x) direction
    Rinv.x = solve(chol(crossprod(B.x) + 1e-10))    
    svdRPR.x = svd(crossprod(Rinv.x, P.x %*% Rinv.x))
    RinvU.x = Rinv.x %*% svdRPR.x$u
    tau = svdRPR.x$d
    if (q.x > 0) tau[(k.x-q.x+1) : k.x] = 0
    A.x = B.x %*% RinvU.x
    
    list(basis.f=basis.f, basis.x=basis.x, B.f=B.f, P.f=P.f, B.x=B.x, P.x=P.x, 
         RinvU.x=RinvU.x, A.x=A.x, tau=tau)
}

# Ledoit-Wolf test
lw.test  = function(X) 
{
  n = nrow(X) - 1
  p = ncol(X)
  S = cov(X)
  trS = sum(diag(S))
  trS2 = sum(S^2)
  stat = (n * p * trS2/trS^2 - n - p - 1)/2
  pvalue = 2 * pnorm(-abs(stat))
  list(stat = stat, pvalue = pvalue)
}

# Newton-Cotes weights
ncw = function(degree, d=10000) {
	x. = seq(0,1,,degree+1)
	w = c()
	for (j in 1:(degree+1)) {
		lj <- function(x) prod((rep(x,degree)-x.[-j])/(rep(x.[j],degree)-x.[-j]))
		ljvec = Vectorize(lj)
		xx = seq(0,1,,d+1) - .5/d
		w[j] = sum(ljvec(xx[-1])) / d
	}
	w
}

# E.g....
# ncw(8)*28350

#' @importFrom fda eval.basis 
pen.nc = function(basis, degree=8, pen.order=2) {
	weights = switch(as.character(degree),
	    "2" = c(1,4,1) / 6,
	    "4" = c(7,32,12,32,7) / 90,
	    "6" = c(41, 216, 27, 272, 27, 216, 41) / 840,
	    "8" = c(989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989) / 28350
	)
	nbasis = basis$nbasis
	lox = c(basis$rangeval[1], basis$param, basis$rangeval[2])
	nlox = length(lox)
	pts = basis$rangeval[1]
	for (k in 1:(nlox-1)) pts = c(pts, seq(lox[k],lox[k+1],,degree+1)[-1])
	f.pts = eval.basis(pts, basis, pen.order)
	penmat = matrix(0, nbasis, nbasis)
	# The following is NOT efficient!
	for (i in 1:nbasis) for (j in 1:nbasis) for (k in 1:(nlox-1)) {		
		fvals = f.pts[(degree*(k-1)+1):(degree*k+1), i] * f.pts[(degree*(k-1)+1):(degree*k+1), j]
		penmat[i,j] = penmat[i,j] + crossprod(weights, fvals) * (lox[k+1]-lox[k])  
	}
	penmat
}

# Knots of smallbasis should be subset of those of bigbasis
#' @importFrom fda eval.basis 
multipen.nc = function(bigbasis, smallbasis, degree=8, pen.order=0) {
	weights = switch(as.character(degree),
	    "2" = c(1,4,1) / 6,
	    "4" = c(7,32,12,32,7) / 90,
	    "6" = c(41, 216, 27, 272, 27, 216, 41) / 840,
	    "8" = c(989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989) / 28350
	)
	nbig = bigbasis$nbasis; nsmall = smallbasis$nbasis
	lox = c(bigbasis$rangeval[1], bigbasis$param, bigbasis$rangeval[2])
	nlox = length(lox)
	pts = bigbasis$rangeval[1]
	for (k in 1:(nlox-1)) pts = c(pts, seq(lox[k],lox[k+1],,degree+1)[-1])
	bigvals = eval.basis(pts, bigbasis, pen.order)
	smallvals = eval.basis(pts, smallbasis, pen.order)
	penmats = vector("list", nsmall)
	for (m in 1:nsmall) {
	    penmats[[m]] = matrix(0, nbig, nbig)
	    for (i in 1:nbig) for (j in 1:nbig) for (k in 1:(nlox-1)) {		
		    fvals = smallvals[(degree*(k-1)+1):(degree*k+1), m] * bigvals[(degree*(k-1)+1):(degree*k+1), i] * bigvals[(degree*(k-1)+1):(degree*k+1), j]
		    penmats[[m]][i,j] = penmats[[m]][i,j] + crossprod(weights, fvals) * (lox[k+1]-lox[k]) 
		}
	}
	penmats
}
