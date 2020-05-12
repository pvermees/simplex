#' @rdname calplot
#' @export
calplot <- function(x,...){ UseMethod("calplot",x) }
#' @rdname calplot
#' @export
calplot.default <- function(x,...){ stop('No default method.') }
#' @rdname calplot
#' @export
calplot.simplex <- function(x,a,b,X=c('UO2','U238'),
                            Y=c('Pb206','U238'),raw=FALSE,...){
    B <- beta2york(b=b,X=X,Y=Y)
    if (raw){
        XX <- NULL
        YY <- NULL
        snames <- names(x)
        for (sname in snames){
            spot <- x[[sname]]
            NXp <- betapars(spot=spot,ion=X[1],a=a[[sname]])
            DXp <- betapars(spot=spot,ion=X[2],a=a[[sname]])
            NYp <- betapars(spot=spot,ion=Y[1],a=a[[sname]])
            DYp <- betapars(spot=spot,ion=Y[2],a=a[[sname]])
            XX <- cbind(XX,log(NXp$sig - NXp$bkg) - log(DXp$sig - DXp$bkg))
            YY <- cbind(YY,log(NYp$sig - NYp$bkg) - log(DYp$sig - DYp$bkg))
        }
        Xlim <- rep(0,2)
        Ylim <- rep(0,2)
        Xlim[1] <- min(B[,'X']-2*B[,'sX'],XX)
        Xlim[2] <- max(B[,'X']+2*B[,'sX'],XX)
        Ylim[1] <- min(B[,'Y']-2*B[,'sY'],YY)
        Ylim[2] <- max(B[,'Y']+2*B[,'sY'],YY)
        IsoplotR:::scatterplot(B,xlim=Xlim,ylim=Ylim)
        matlines(XX,YY,lty=1,col='grey')
        #points(XX[1,],YY[1,],pch=21,bg='black')
        #points(XX[nrow(XX),],YY[nrow(YY),],pch=21,bg='white')
    } else {
        IsoplotR:::scatterplot(B,...)
    }
}
#' @rdname calplot
#' @export
calplot.standards <- function(x,a,b,X=c('UO2','U238'),
                              Y=c('Pb206','U238'),raw=FALSE,...){
    calplot(x=x$x,a=a,b=b,X=X,Y=Y,raw=raw,...)
}

beta2york <- function(b,X=c('UO2','U238'),Y=c('Pb206','U238')){
    snames <- names(b)
    ns <- length(snames)
    out <- matrix(0,nrow=ns,ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        b0g <- b[[sname]]$b0g
        covmat <- b[[sname]]$cov
        Xlab <- paste0('b0[',X[1],'/',X[2],']')
        Ylab <- paste0('b0[',Y[1],'/',Y[2],']')
        out[sname,'X'] <- b0g[Xlab]
        out[sname,'Y'] <- b0g[Ylab]
        out[sname,'sX'] <- sqrt(covmat[Xlab,Xlab])
        out[sname,'sY'] <- sqrt(covmat[Ylab,Ylab])
        out[sname,'rXY'] <- covmat[Xlab,Ylab]/(out[sname,'sX']*out[sname,'sY'])
    }
    out
}
