calibration <- function(lr,dc=NULL,dat=NULL,type="U-Pb",oxide=NULL,plot=1,t=0,...){
    if (stable(type)){
        fit <- stable_calibration(lr=lr,dc=dc,dat=dat,type=type,plot=plot,...)
    } else {
        fit <- geochron_calibration(lr=lr,dc=dc,dat=dat,type=type,
                                    oxide=oxide,plot=plot,t=t,...)
    }
    invisible(fit)
}

stable_calibration <- function(lr,dc=NULL,dat=NULL,
                               type="dO18",plot=1,...){
    LL <- function(par,lr=lr){
        out <- 0
        snames <- names(lr)
        for (sname in snames){
            LL <- mahalanobis(x=lr[[sname]],center=init,cov=lr[[sname]]$cov)
            out <- out + LL
        }
        -out
    }
    init <- lr[[1]]$lr
    wtdmean <- optim(init,fn=LL,gr=NULL,method='BFGS',hessian=TRUE,lr=lr)
    out <- list()
    out$x <- rep(0,nn)
    out$cov <- matrix(0,nn,nn)
    out
}

geochron_calibration <- function(lr,dc=NULL,dat=NULL,type="U-Pb",
                                 oxide=NULL,plot=1,t=0,...){
    if (type=="U-Pb"){
        if (is.null(oxide)){
            b0gnames <- names(lr[[1]]$b0g)
            if (any(grepl('UO2',b0gnames))){
                oxide <- 'UO2'
            } else if (any(grepl('UO',b0gnames))){
                oxide <- 'UO'
            } else {
                stop('No valid oxide was measured.')
            }
        }
        x <- c(oxide,'U238')
        y <- c('Pb206','U238')
    }
    B <- beta2york(lr=lr,t=t,x=x,y=y)
    if (plot>0){
        xlab <- paste0('log[',x[1],'/',x[2],']')
        ylab <- paste0('log[',y[1],'/',y[2],']')
        if (plot==1){
            fit <- IsoplotR:::isochron(B,xlab=xlab,ylab=ylab,...)
        } else if (!is.null(dat)){
            X <- NULL
            Y <- NULL
            snames <- names(dat)
            for (sname in snames){
                sp <- spot(dat=dat,sname=sname)
                Nxp <- betapars(spot=sp,ion=x[1],dc=dc[[sname]])
                Dxp <- betapars(spot=sp,ion=x[2],dc=dc[[sname]])
                Nyp <- betapars(spot=sp,ion=y[1],dc=dc[[sname]])
                Dyp <- betapars(spot=sp,ion=y[2],dc=dc[[sname]])
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
        } else {
            stop("Can't plot the data.")
        }
    } else {
        fit <- IsoplotR:::york(B)
    }
    fit
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
