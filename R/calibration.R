calibration <- function(lr,dc=NULL,dat=NULL,oxide=NULL,plot=1,t=0,...){
    if (stable(lr)){
        fit <- stable_calibration(lr=lr,dc=dc,dat=dat,type=type,plot=plot,...)
    } else {
        fit <- geochron_calibration(lr=lr,dc=dc,dat=dat,
                                    oxide=oxide,plot=plot,t=t,...)
    }
    invisible(fit)
}

stable_calibration <- function(lr,dc=NULL,dat=NULL,plot=1,...){
    LL <- function(par,LR){
        out <- 0
        np <- length(par)
        snames <- names(LR)
        for (sname in snames){
            X <- LR[[sname]]$b0g[1:np]
            E <- LR[[sname]]$cov[1:np,1:np]
            LL <- mahalanobis(x=X,center=par,cov=E)
            out <- out + LL
        }
        out
    }
    ni <- length(lr$num)
    LR <- lr$x
    init <- LR[[1]]$b0g[1:ni]
    wtdmean <- optim(init,fn=LL,gr=NULL,method='BFGS',hessian=TRUE,LR=LR)
    out <- list()
    out$x <- wtdmean$par
    out$cov <- solve(wtdmean$hessian)
    if (plot>0){
        np <- length(lr$num)-1    # number of plot panels
        nr <- ceiling(sqrt(np)) # number of rows
        nc <- ceiling(np/nr)    # number of columns
        oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
        for (i in 1:nr){
            for (j in (i+1):max(nc,nr+1)){
                b0names <- names(out$x)
                xlab <- paste0(lr$num[i],'/',lr$den[i])
                ylab <- paste0(lr$num[j],'/',lr$den[j])
                B <- beta2york(lr=lr,
                               x=c(lr$num[i],lr$den[i]),
                               y=c(lr$num[j],lr$den[j]))
                IsoplotR::scatterplot(B)
                graphics::mtext(side=1,text=xlab,line=2)
                graphics::mtext(side=2,text=ylab,line=2)
                ell <- IsoplotR::ellipse(out$x[i],out$x[j],out$cov[c(i,j),c(i,j)])
                graphics::polygon(ell,col='white')
            }
        }
    }
    out
}

geochron_calibration <- function(lr,dc=NULL,dat=NULL,oxide=NULL,plot=1,t=0,...){
    if (datatype(lr)=="U-Pb"){
        if (is.null(oxide)){
            if ('UO2'%in%lr$num){
                oxide <- 'UO2'
            } else if ('UO'%in%lr$num){
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
        b <- b0gt2b(LR=lr$x[[sname]],t=t,
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

b0gt2b <- function(LR,t=0,xlab,ylab){
    b0g <- LR$b0g
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
    out$cov <- J %*% LR$cov %*% t(J)
    out
}
