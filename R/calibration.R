calibration <- function(lr,dc=NULL,dat=NULL,x=c('UO2','U238'),
                        y=c('Pb206','U238'),plot=1,t=0,...){
    B <- beta2york(lr=lr,t=t,x=x,y=y)
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

beta2york <- function(lr,t=0,x=c('UO2','U238'),y=c('Pb206','U238')){
    snames <- names(lr)
    ns <- length(snames)
    out <- matrix(0,nrow=ns,ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        b <- b0gt2b(lr=lr[[sname]],t=t,
                    xlab=paste0(x[1],'/',x[2]),
                    ylab=paste0(y[1],'/',y[2]))
        out[sname,'X'] <- b$x[1]
        out[sname,'Y'] <- b$x[2]
        out[sname,'sX'] <- sqrt(b$cov[1,1])
        out[sname,'sY'] <- sqrt(b$cov[2,2])
        out[sname,'rXY'] <- b$cov[1,2]/(out[sname,'sX']*out[sname,'sY'])
    }
    out
}

b0gt2b <- function(lr,t=0,xlab,ylab){
    b0g <- lr$b0g
    b0gnames <- names(b0g)
    nb0g <- length(b0gnames)
    out <- list()
    out$x <- rep(0,2)
    J <- matrix(0,2,nb0g)
    ibx <- which(b0gnames %in% paste0('b0[',xlab,']'))
    iby <- which(b0gnames %in% paste0('b0[',ylab,']'))
    igx <- which(b0gnames %in% paste0('g[',xlab,']'))
    igy <- which(b0gnames %in% paste0('g[',ylab,']'))
    out$x[1] <- b0g[ibx] + b0g[igx]*hours(t)
    out$x[2] <- b0g[iby] + b0g[igy]*hours(t)
    J[1,ibx] <- 1
    J[1,igx] <- hours(t)
    J[2,iby] <- 1
    J[2,igy] <- hours(t)
    out$cov <- J %*% lr$cov %*% t(J)
    out
}
