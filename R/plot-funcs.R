#' @title Plot a vsm2s object
#' 
#' @description Plot function for objects created by \code{\link{vsm2s}}
#' 
#' @param x object created by \code{\link{vsm2s}}
#' @param newx new predictor values (at which to predict)
#' @param xlab,ylab labels for x-axis and y-axis
#' @param label.curves whether to label the curves
#' @param argvals function argument values 
#' @param which which of the fits to plot (relevant if \code{x} includes multiple fits)
#' @param lty line type
#' @param owncolor user-define color scheme
#' @param ... other arguments, passed to \code{\link[graphics]{matplot}}
#' @importFrom colorRamps blue2red
#' @importFrom grDevices colorRamp rgb
#' 
#' @export

plot.vsm2s = function(x, newx=sort(x$x), xlab="", ylab="", 
                      label.curves=FALSE, argvals=x$argvals, 
                      which=1, lty=1, owncolor=colorRamp(blue2red(100)), ...) { 
    nx = length(newx)
	predv = predict(x, newx, which=which)
	colvec = rgb(owncolor((newx-newx[1])/(newx[nx]-newx[1])),maxColorValue=255)
	matplot(x$argvals, t(predv), col=colvec, type="l", lty=lty, xlab=xlab, ylab=ylab, ...)
    if (label.curves) text(rep(quantile(argvals, .95), length(newx)), predv[,which.min(abs(argvals - quantile(argvals, .95)))], newx, col=colvec)
}

#' @title Plot a tps.fd object
#' 
#' @description Plot function for objects created by \code{\link{tps.fd}}
#' 
#' @param x object created by \code{\link{tps.fd}}
#' @param newx new predictor values (at which to predict)
#' @param xlab,ylab labels for x-axis and y-axis
#' @param label.curves whether to label the curves
#' @param argvals function argument values 
#' @param lty line type
#' @param owncolor user-defined color scheme
#' @param ... other arguments, passed to \code{\link[graphics]{matplot}}
#' @importFrom colorRamps blue2red
#' @importFrom grDevices colorRamp rgb
#' 
#' @export

plot.tps.fd = function(x, newx=sort(x$x), xlab="", ylab="", label.curves=FALSE, argvals=x$argvals, lty=1, owncolor=colorRamp(blue2red(100)), ...) { 
    nx = length(newx)
    predv = predict(x, newx)
    colvec = rgb(owncolor((newx-newx[1])/(newx[nx]-newx[1])),maxColorValue=255)
    matplot(x$argvals, t(predv), col=colvec, type="l", lty=lty, xlab=xlab, ylab=ylab, ...)
    if (label.curves) text(rep(quantile(argvals, .95), length(newx)), predv[,which.min(abs(argvals - quantile(argvals, .95)))], newx, col=colvec)
}

#' @title Rainbow plot  
#' 
#' @description Rainbow plot for functional data, color-coded by a predictor
#' 
#' @param Y n x L matrix of functional responses
#' @param x n-dimensional predictor vector 
#' @param argvals function argument values
#' @param subset subset of functional observations to plot
#' @param label.curves whether to label the curves
#' @param xlab,ylab labels for x-axis and y-axis
#' @param legend logical: Should a legend be added?
#' @param digit,posx,posy,zval passed to \code{\link[shape]{colorlegend}}
#' @param ... other arguments, passed to \code{\link{matplot}}
#' @importFrom grDevices colorRamp rgb
#' @importFrom colorRamps blue2red
#' 
#' @export

rainbow.plot = function(Y, x, argvals=NULL, subset = 1:length(x), label.curves=FALSE, xlab="", ylab="", legend=TRUE, digit=1, posx=c(.9,.93), posy=c(0.05,0.9), zval=NULL,...)   {
    nx = length(subset)
    newx = x[subset] 
    predv = Y[subset, ] 
    if (is.null(argvals)) argvals = 1:dim(predv)[2]
    owncolor = colorRamp(blue2red(100))
    colvec = rgb(owncolor((newx-max(newx))/(min(newx)-max(newx))),maxColorValue=255)
    
    matplot(argvals, t(predv), col=colvec, type="l", lty=1, xlab=xlab, ylab=ylab, ...)
    if (label.curves) text(rep(quantile(argvals, .95), length(newx)), predv[,which.min(abs(argvals - quantile(argvals, .95)))], newx, col=colvec)
    if (legend) {
    	    legvec = rgb(owncolor(sort((newx-max(newx))/(min(newx)-max(newx)), TRUE)),maxColorValue=255)
    	    shape::colorlegend(col=legvec, zlim=range(newx), main="", sub="x", left=FALSE, digit=digit, posx=posx, posy=posy, zval=zval)
    }
}
