calibration <- function(lr,stand,prefix=NULL,snames=NULL,i=NULL,
                        invert=FALSE,oxide=NULL,t=0){
    out <- lr
    dat <- subset(x=lr,prefix=prefix,snames=snames,i=i,invert=invert)
    if (stable(lr)) out$cal <- stable_calibration(lr=dat)
    else out$cal <- geochron_calibration(lr=dat,oxide=oxide,t=t)
    out$stand <- stand
    class(out) <- append('calibration',class(out))
    out
}

stable_calibration <- function(lr){
    LL <- function(par,lr){
        out <- 0
        np <- length(par)
        snames <- names(lr$x)
        for (sname in snames){
            X <- lr$x[[sname]]$lr$b0g[1:np]
            E <- lr$x[[sname]]$lr$cov[1:np,1:np]
            LL <- mahalanobis(x=X,center=par,cov=E)
            out <- out + LL
        }
        out
    }
    ni <- length(lr$num)
    init <- lr$x[[1]]$lr$b0g[1:ni]
    wtdmean <- optim(init,fn=LL,gr=NULL,method='BFGS',hessian=TRUE,lr=lr)
    out <- list()
    out$x <- wtdmean$par
    out$cov <- solve(wtdmean$hessian)
    out
}

geochron_calibration <- function(lr,oxide=NULL,t=0,...){
    out <- list()
    out$t <- t
    type <- datatype(lr)
    if (type=="U-Pb"){
        out[[type]] <- add_geochron_calibration(lr,oxide=oxide,t=t,type)
    } else if (datatype(lr)=="U-Th-Pb"){
        out[['U-Pb']] <- add_geochron_calibration(lr,oxide=oxide,t=t,type="U-Pb")
        out[['Th-Pb']] <- add_geochron_calibration(lr,oxide=oxide,t=t,type="Th-Pb")
    } else {
        stop("Invalid data type.")
    }
    out
}

add_geochron_calibration <- function(lr,oxide=NULL,t=0,type='U-Pb'){
    if (type=='U-Pb'){
        if (is.null(oxide)){
            if ('UO2'%in%lr$num){
                oxide <- 'UO2'
            } else if ('UO'%in%lr$num){
                oxide <- 'UO'
            } else {
                stop('No valid oxide was measured.')
            }
        }
        num <- c(oxide,'Pb206')
        den <- c('U238','U238')
    } else if (type=='Th-Pb'){
        if ('ThO2'%in%lr$num){
            oxide <- 'ThO2'
        } else if ('ThO'%in%lr$num){
            oxide <- 'ThO'
        } else {
            stop('No valid Th oxide was measured.')
        }
        num <- c(oxide,'Pb208')
        den <- c('Th232','Th232')
    } else {
        stop("Invalid data type.")
    }
    out <- list(num=num,den=den,oxide=oxide)
    out$york <- beta2york(lr=lr,t=t,num=num,den=den)
    out$fit <- IsoplotR:::york(out$york)
    out
}

beta2york <- function(lr,t=0,num=c('UO2','Pb206'),
                      den=c('U238','U238'),snames=NULL,i=NULL){
    if (is.null(snames)) snames <- names(lr$x)
    if (!is.null(i)) snames <- snames[i]
    ns <- length(snames)
    out <- matrix(0,nrow=ns,ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        b <- b0gt2b(LR=lr$x[[sname]]$lr,t=t,
                    xlab=paste0(num[1],'/',den[1]),
                    ylab=paste0(num[2],'/',den[2]))
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

plot.calibration <- function(cal,option=1,snames=NULL,i=NULL,...){
    if (stable(cal)) calplot_stable(dat=cal,snames=snames,i=i,...)
    else calplot_geochronology(dat=cal,option=option,snames=snames,i=i,...)
}

calplot_stable <- function(dat,snames=NULL,i=NULL,...){
    np <- length(dat$num)-1  # number of plot panels
    nr <- ceiling(sqrt(np))  # number of rows
    nc <- ceiling(np/nr)     # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (i in 1:nr){
        for (j in (i+1):max(nc,nr+1)){
            b0names <- names(dat$x)
            xlab <- paste0(dat$num[i],'/',dat$den[i])
            ylab <- paste0(dat$num[j],'/',dat$den[j])
            B <- beta2york(lr=dat,snames=snames,i=i,
                           num=dat$num[c(i,j)],den=dat$den[c(i,j)])
            IsoplotR::scatterplot(B,...)
            graphics::mtext(side=1,text=xlab,line=2)
            graphics::mtext(side=2,text=ylab,line=2)
            ell <- IsoplotR::ellipse(dat$cal$x[i],dat$cal$x[j],
                                     dat$cal$cov[c(i,j),c(i,j)])
            graphics::polygon(ell,col='white')
        }
    }
}

calplot_geochronology <- function(dat,option=1,snames=NULL,i=NULL,...){
    if (is.null(snames)) snames <- rownames(dat$cal$york)
    if (!is.null(i)) snames <- rownames(dat$cal$york)[i]
    yd <- dat$cal$york[snames,,drop=FALSE]
    num <- dat$cal$num
    den <- dat$cal$den
    xlab <- paste0('log[',num[1],'/',den[1],']')
    ylab <- paste0('log[',num[2],'/',den[2],']')
    if (option==1){
        IsoplotR:::scatterplot(yd,fit=dat$cal$fit,...)
    } else {
        X <- NULL
        Y <- NULL
        for (sname in snames){
            sp <- spot(dat=dat,sname=sname)
            Nxp <- betapars(spot=sp,ion=num[1])
            Dxp <- betapars(spot=sp,ion=den[1])
            Nyp <- betapars(spot=sp,ion=num[2])
            Dyp <- betapars(spot=sp,ion=den[2])
            newX <- Nxp$g*(Dxp$t-Nxp$t) +
                log(Nxp$sig - Nxp$bkg) - log(Dxp$sig - Dxp$bkg)
            newY <- Nyp$g*(Dyp$t-Nyp$t) +
                log(Nyp$sig - Nyp$bkg) - log(Dyp$sig - Dyp$bkg)
            X <- cbind(X,newX)
            Y <- cbind(Y,newY)
        }
        Xlim <- rep(0,2)
        Ylim <- rep(0,2)
        Xlim[1] <- min(yd[,'X']-2*yd[,'sX'],X)
        Xlim[2] <- max(yd[,'X']+2*yd[,'sX'],X)
        Ylim[1] <- min(yd[,'Y']-2*yd[,'sY'],Y)
        Ylim[2] <- max(yd[,'Y']+2*yd[,'sY'],Y)
        IsoplotR:::scatterplot(yd,fit=dat$cal$fit,xlim=Xlim,
                               ylim=Ylim,xlab=xlab,ylab=ylab,...)
        matlines(X,Y,lty=1,col='darkgrey')
        if (option>2){
            points(X[1,],Y[1,],pch=21,bg='black')
            points(X[nrow(X),],Y[nrow(Y),],pch=21,bg='white')
        }
    }
}
