calibration <- function(lr,dc,dat=NULL,x=c('UO2','U238'),
                                y=c('Pb206','U238'),plot=1,...){
    B <- beta2york(lr=lr,x=x,y=y)
    fit <- IsoplotR:::york(B)
    if (plot>0){
        xlab <- paste0('log[',x[1],'/',x[2],']')
        ylab <- paste0('log[',y[1],'/',y[2],']')
    }
    if (plot==1){
        IsoplotR:::isochron(B,xlab=xlab,ylab=ylab,...)
    }
    if (plot>1 & !is.null(dat)){
        if (methods::is(dat,'standards')) d <- dat$x
        else d <- dat
        X <- NULL
        Y <- NULL
        snames <- names(d)
        for (sname in snames){
            spot <- d[[sname]]
            Nxp <- betapars(spot=spot,ion=x[1],dc=dc[[sname]])
            Dxp <- betapars(spot=spot,ion=x[2],dc=dc[[sname]])
            Nyp <- betapars(spot=spot,ion=y[1],dc=dc[[sname]])
            Dyp <- betapars(spot=spot,ion=y[2],dc=dc[[sname]])
            X <- cbind(X,log(Nxp$sig - Nxp$bkg) - log(Dxp$sig - Dxp$bkg))
            Y <- cbind(Y,log(Nyp$sig - Nyp$bkg) - log(Dyp$sig - Dyp$bkg))
        }
        Xlim <- rep(0,2)
        Ylim <- rep(0,2)
        Xlim[1] <- min(B[,'X']-2*B[,'sX'],X)
        Xlim[2] <- max(B[,'X']+2*B[,'sX'],X)
        Ylim[1] <- min(B[,'Y']-2*B[,'sY'],Y)
        Ylim[2] <- max(B[,'Y']+2*B[,'sY'],Y)
        fit <- IsoplotR:::isochron(B,xlim=Xlim,ylim=Ylim,xlab=xlab,ylab=ylab)
        matlines(X,Y,lty=1,col='darkgrey')
        if (plot>2){
            points(X[1,],Y[1,],pch=21,bg='black')
            points(X[nrow(X),],Y[nrow(Y),],pch=21,bg='white')
        }
    }
    invisible(fit)
}

beta2york <- function(lr,x=c('UO2','U238'),y=c('Pb206','U238')){
    snames <- names(lr)
    ns <- length(snames)
    out <- matrix(0,nrow=ns,ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        b0g <- lr[[sname]]$b0g
        covmat <- lr[[sname]]$cov
        Xlab <- paste0('b0[',x[1],'/',x[2],']')
        Ylab <- paste0('b0[',y[1],'/',y[2],']')
        out[sname,'X'] <- b0g[Xlab]
        out[sname,'Y'] <- b0g[Ylab]
        out[sname,'sX'] <- sqrt(covmat[Xlab,Xlab])
        out[sname,'sY'] <- sqrt(covmat[Ylab,Ylab])
        out[sname,'rXY'] <- covmat[Xlab,Ylab]/(out[sname,'sX']*out[sname,'sY'])
    }
    out
}
