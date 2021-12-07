#' @title calibrate SIMS data
#' @description convert signal logratios to atomic logratios
#' @param cal an object of class \code{calibration}
#' @param exterr logical flag indicating whether the uncertainty of
#'     the calibration should be propagated.
#' @return an object of class \code{calibrated}
#' @examples
#' data('Cameca_oxygen',package='simplex')
#' dc <- drift(x=Cameca_oxygen)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28'))
#' cd <- calibrate(cal)
#' del <- delta(cd,log=FALSE)
#' tab <- data2table(del)
#' @export
calibrate <- function(cal,exterr=FALSE){
    if (ncol(cal$calibration$pairing)==2){
        out <- calibrate_average(dat=cal,exterr=exterr)
    } else {
        out <- calibrate_regression(dat=cal,exterr=exterr)
    }
    out
}

calibrate_average <- function(dat,exterr=FALSE){
    out <- dat
    cal <- dat$calibration
    ns <- length(dat$samples)
    tavg <- time_average(dat,t=cal$t)
    alliso <- names(tavg[[1]]$val)
    caliso <- cal$cal$ratios
    staiso <- cal$stand$ratios
    nai <- length(alliso)
    nci <- length(caliso)
    nsi <- length(staiso)
    E <- matrix(0,nrow=ns*nai+nsi+nci,ncol=ns*nai+nsi+nci)
    J <- matrix(0,nrow=ns*nai,ncol=ns*nai+nsi+nci)
    i1 <- lookup(alliso,cal$pairing$smp)
    i2 <- lookup(caliso,alliso)
    E[ns*nai+1:nsi,ns*nai+1:nsi] <- as.matrix(cal$stand[,-c(1,2)])
    E[ns*nai+nsi+1:nci,ns*nai+nsi+1:nci] <- as.matrix(cal$cal[,-c(1,2)])
    val <- rep(0,ns*nai)
    for (i in 1:ns){
        i3 <- (i-1)*nai+1:nai
        val[i3] <- tavg[[i]]$val
        if (exterr){
            val[i3[i2]] <- val[i3[i2]] + cal$stand$val[i1] - cal$cal$val
            J[i3[i1],ns*nai+1:nsi] <- diag(nci)
            J[i3[i2],ns*nai+nsi+1:nci] <- -diag(nci)
        }
        E[i3,i3] <- tavg[[i]]$cov
        J[i3,i3] <- diag(nai)
    }
    out$calibrated <- list(snames=names(dat$samples),
                           ratios=paste0(dat$method$num,'/',dat$method$den),
                           val=val,cov=J %*% E %*% t(J))
    class(out) <- unique(append("calibrated",class(out)))
    out
}

calibrate_regression <- function(dat,exterr=FALSE){
    out <- dat
    cal <- dat$calibration
    ns <- length(dat$samples)
    tavg <- time_average(dat,t=cal$t)
    ncal <- nrow(cal$cal)
    alliso <- names(tavg[[1]]$val)
    outiso <- alliso[!(alliso %in% cal$pairing$versus)]
    staiso <- cal$stand$ratios
    nai <- length(alliso)
    noi <- length(outiso)
    nsi <- length(staiso)
    nab <- nrow(cal$pairing)
    iout <- lookup(outiso,alliso)
    istd <- ns*nai+1:nsi
    iab <- matrix(0,nab,2)
    iDPall <- rep(NA,nab)
    iDPout <- rep(NA,nab)
    iDPstd <- rep(NA,nab)
    iOP <- rep(NA,nab)
    for (i in 1:nab){
        iab[i,] <- ns*nai+nsi+(i-1)*2+1:2
        iDPall[i] <- lookup(cal$cal[i,'y'],alliso)
        iDPout[i] <- lookup(cal$cal[i,'y'],outiso)
        iDPstd[i] <- lookup(cal$cal[i,'y'],staiso)
        iOP[i] <- lookup(cal$cal[i,'x'],alliso)
    }
    val <- rep(0,ns*noi)
    E <- matrix(0,nrow=ns*nai+nsi+2*nab,ncol=ns*nai+nsi+2*nab)
    E[istd,istd] <- as.matrix(cal$stand[1:nsi,2+1:nsi])
    for (i in 1:nab){
        E[iab[i,1],iab[i,1]] <- cal$cal[i,'s[a]']^2
        E[iab[i,2],iab[i,2]] <- cal$cal[i,'s[b]']^2
        E[iab[i,1],iab[i,2]] <- cal$cal[i,'cov.ab']
        E[iab[i,2],iab[i,1]] <- cal$cal[i,'cov.ab']
    }
    J <- matrix(0,nrow=ns*noi,ncol=ns*nai+nsi+2*nab)
    for (i in 1:ns){
        ioi <- (i-1)*noi+1:noi
        val[ioi] <- tavg[[i]]$val[iout]
        iai <- (i-1)*nai+1:nai
        E[iai,iai] <- tavg[[i]]$cov
        J[ioi,iai] <- diag(nai)[iout,] # dval/dval
        for (j in 1:nab){
            val[ioi[iDPout[j]]] <- tavg[[i]]$val[iDPall[j]] -
                cal$cal[j,'a'] - cal$cal[j,'b']*tavg[[i]]$val[iOP[j]] +
                cal$stand[iDPstd[j],'val']
            J[ioi[iDPout[j]],iai[iOP[j]]] <- - cal$cal[j,'b'] # dval/dOP
            if (exterr){
                J[ioi[iDPout[j]],iab[j,1]] <- -1 # dval/da
                J[ioi[iDPout[j]],iab[j,2]] <- -tavg[[i]]$val[iOP[j]] # dval/db
                J[ioi[iDPout[j]],istd[iDPstd[j]]] <- 1 # dval/dDPstd
            }
        }
    }
    out$calibrated <- list(snames=names(dat$samples),ratios=outiso,
                           val = val, cov = J%*% E %*% t(J))
    out
}

#' @title plot calibrated data
#' @description shows the calibrated data on a logratio plot.
#' @param x an object of class \code{calibrated}
#' @param option for U-Pb or Th-Pb data. If \code{1}, plots the data
#'     as error ellipses; if \code{2}, adds the raw data; if \code{3},
#'     marks the first and last measurement with a black and white
#'     circle, respectively.
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function
#' @examples
#' data('Cameca_UPb')
#' dc <- drift(x=Cameca_UPb)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='Plesovice'))
#' cd <- calibrate(cal)
#' plot(cd)
#' @method plot calibrated
#' @export
plot.calibrated <- function(x,option=1,...){
    if (ncol(x$calibration$pairing)==2){
        out <- caldplot_stable(dat=x,...)
    } else {
        out <- caldplot_geochronology(dat=x,option=option,...)
    }
    invisible(out)
}

caldplot_stable <- function(dat,...){
    snames <- names(dat$samples)
    tab <- data2table.calibrated(dat)
    nn <- length(dat$calibrated$ratios)
    if (nn>1){
        np <- nn*(nn-1)/2       # number of plot panels
        nc <- ceiling(sqrt(np)) # number of rows
        nr <- ceiling(np/nc)    # number of columns
        oldpar <- graphics::par(mfrow=c(nr,nc),mar=rep(3.5,4))
        ii <- 1
        for (i in 1:(nn-1)){
            for (j in (i+1):nn){
                xlab <- dat$calibrated$ratios[i]
                ylab <- dat$calibrated$ratios[j]
                X <- tab[,xlab]
                Y <- tab[,ylab]
                sX <- tab[,paste0('s[',xlab,']')]
                sY <- tab[,paste0('s[',ylab,']')]
                rXY <- tab[,paste0('r[',xlab,',',ylab,']')]
                calval <- dat$calibration$cal
                calcov <- dat$calibration$cal[,-c(1:2)]
                ical <- match(xlab,calval$ratios)
                jcal <- match(ylab,calval$ratios)
                x <- calval[ical,'val']
                sx <- sqrt(calcov[ical,ical])
                y <- calval[jcal,'val']
                sy <- sqrt(calcov[jcal,jcal])
                xlim <- c(min(x-3*sx,X-3*sX),max(x+3*sx,X+3*sX))
                ylim <- c(min(y-3*sy,Y-3*sY),max(y+3*sy,Y+3*sY))
                graphics::plot(xlim,ylim,type='n',ann=FALSE)
                deltagrid(dat,ical,jcal)
                IsoplotR::scatterplot(cbind(X,sX,Y,sY,rXY),
                                      xlim=xlim,ylim=ylim,add=TRUE,...)
                graphics::mtext(side=1,text=xlab,line=2)
                graphics::mtext(side=2,text=ylab,line=2)
                ell <- IsoplotR::ellipse(x,y,calcov[c(ical,jcal),c(ical,jcal)])
                graphics::polygon(ell,col='white')
                ii <- ii + 1
                if (ii>np) break
            }
        }
        graphics::par(oldpar)
    } else {
        ns <- length(snames)
        tfact <- qnorm(0.975)
        lr <- tab[,1]
        ll <- lr - tfact*tab[,2]
        ul <- lr + tfact*tab[,2]
        matplot(rbind(1:ns,1:ns),rbind(ll,ul),type='l',lty=1,
                col='black',bty='n',xlab='sample #',ylab='')
        points(1:ns,lr,pch=16)
        cal <- dat$calibration$cal
        lines(c(1,ns),rep(cal$val,2),lty=2)
        lines(c(1,ns),rep(cal$val-tfact*sqrt(cal$cov),2),lty=3)
        lines(c(1,ns),rep(cal$val+tfact*sqrt(cal$cov),2),lty=3)
    }
}

deltagrid <- function(dat,ical,jcal){
    usr <- graphics::par('usr')
    xlim <- usr[1:2]
    ylim <- usr[3:4]
    ratios <- dat$calibration$cal$ratios
    xs <- dat$calibration$cal$val[ical]
    ys <- dat$calibration$cal$val[jcal]
    ist <- match(ratios[ical],dat$calibration$pairing$smp)
    jst <- match(ratios[jcal],dat$calibration$pairing$smp)
    di <- dat$calibration$stand$val[ist]
    dj <- dat$calibration$stand$val[jst]
    dxmin <- 1000*(xlim[1]-xs)+di
    dxmax <- 1000*(xlim[2]-xs)+di
    dymin <- 1000*(ylim[1]-ys)+dj
    dymax <- 1000*(ylim[2]-ys)+dj
    dxticks <- pretty(c(dxmin,dxmax))
    dyticks <- pretty(c(dymin,dymax))
    xticks <- (dxticks-di)/1000+xs
    yticks <- (dyticks-dj)/1000+ys
    nxt <- length(xticks)
    nyt <- length(yticks)
    graphics::matlines(rbind(xticks,xticks),matrix(rep(ylim,nxt),ncol=nxt),
                       lty=3,col='black')
    graphics::matlines(matrix(rep(xlim,nyt),ncol=nyt),
                       rbind(yticks,yticks),lty=3,col='black')
    graphics::axis(side=3,at=xticks,labels=dxticks)
    graphics::mtext(expression(delta*"'"),side=3,line=2)
    graphics::axis(side=4,at=yticks,labels=dyticks)
    graphics::mtext(expression(delta*"'"),side=4,line=2)
    graphics::lines(rep(xs,2),ylim,lty=2,col='red')
    graphics::text(xs,ylim[1],labels=dat$calibration$prefix,pos=4,srt=90,offset=0)
    graphics::lines(xlim,rep(ys,2),lty=2,col='red')
    graphics::text(xlim[1],ys,labels=dat$calibration$prefix,pos=4,offset=0)
}

caldplot_geochronology <- function(dat,option=1,...){
    cal <- dat$calibration
    num <- cal$num
    den <- cal$den
    snames <- names(dat$samples)
    yd <- beta2york(lr=dat,t=seconds(cal$t),snames=snames,num=num,den=den)
    xlab <- paste0('log[',num[1],'/',den[1],']')
    ylab <- paste0('log[',num[2],'/',den[2],']')
    xlim <- rep(0,2)
    xlim[1] <- min(yd[,'X']-3*yd[,'sX'])
    xlim[2] <- max(yd[,'X']+3*yd[,'sX'])
    ylim <- rep(0,2)
    ylim[1] <- min(yd[,'Y']-3*yd[,'sY'],cal$fit$a[1]+cal$fit$b[1]*xlim[1])
    ylim[2] <- max(yd[,'Y']+3*yd[,'sY'],cal$fit$a[1]+cal$fit$b[1]*xlim[2])
    graphics::plot(xlim,ylim,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,type='n')
    agegrid(dat,xlim,ylim)
    if (option==1){
        IsoplotR::scatterplot(yd,fit=cal$fit,add=TRUE,...)
    } else {
        X <- NULL
        Y <- NULL
        for (sname in snames){
            sp <- spot(dat=dat,sname=sname)
            Op <- betapars(spot=sp,ion=num[1])
            Pp <- betapars(spot=sp,ion=den[1])
            Dp <- betapars(spot=sp,ion=num[2])
            Cp <- betapars(spot=sp,ion=num[3])
            CD <- (Cp$sig-Cp$bkg)/(Dp$sig-Dp$bkg)
            CDdc <- exp(Dp$g*(Dp$t-Cp$t))
            newX <- log(Op$sig-Op$bkg) - log(Pp$sig-Pp$bkg) + Op$g*(Pp$t-Op$t)
            newY <- log(Dp$sig-Dp$bkg) - log(Pp$sig-Pp$bkg) + Dp$g*(Pp$t-Dp$t)
            X <- cbind(X,newX)
            Y <- cbind(Y,newY)
        }
        xlim[1] <- min(xlim[1],yd[snames,'X']-3*yd[snames,'sX'],X)
        xlim[2] <- max(xlim[2],yd[snames,'X']+3*yd[snames,'sX'],X)
        ylim[1] <- min(ylim[1],yd[snames,'Y']-3*yd[snames,'sY'],Y)
        ylim[2] <- max(ylim[2],yd[snames,'Y']+3*yd[snames,'sY'],Y)
        IsoplotR::scatterplot(yd,fit=cal$fit,add=TRUE,...)
        graphics::matlines(X,Y,lty=1,col='darkgrey')
        if (option>2){
            graphics::points(X[1,],Y[1,],pch=21,bg='black')
            graphics::points(X[nrow(X),],Y[nrow(Y),],pch=21,bg='white')
        }
    }
}

agegrid <- function(dat,xlim,ylim){
    ratio <- ifelse(datatype(dat)=='U-Pb','Pb206U238','Pb208Th232')
    usr <- graphics::par('usr')
    cal <- dat$calibration
    st <- do.call(dat$standard$fetchfun,args=list(dat=dat))
    lrlim <- rep(0,2)
    a <- cal$fit$a[1]
    b <- cal$fit$b[1]
    adj1 <- st$lr[ratio] - a - b*xlim[2]
    adj2 <- st$lr[ratio] - a - b*xlim[1]
    lrlim[1] <- ylim[1] + adj1
    lrlim[2] <- ylim[2] + adj2
    tlim <- IsoplotR::age(exp(lrlim),method=chronometer(dat))
    tticks <- pretty(tlim)
    nt <- length(tticks)
    lrmin <- log(IsoplotR::age2ratio(tticks,ratio=ratio)[,1]) - adj2
    lrmax <- lrmin + cal$fit$b[1]*diff(xlim)
    xticks <- list(x=NULL,t=NULL)
    yticks <- list(y=NULL,t=NULL)
    for (i in 1:nt){
        if (lrmax[i]>usr[4]){
            xticks$x <- c(xticks$x,(usr[4]-lrmin[i])/b+xlim[1])
            xticks$t <- c(xticks$t,tticks[i])
        } else {
            yticks$y <- c(yticks$y,lrmax[i])
            yticks$t <- c(yticks$t,tticks[i])
        }
    }
    graphics::axis(side=3,at=xticks$x,labels=xticks$t)
    graphics::mtext(side=3,line=2,'t (Ma)')
    graphics::axis(side=4,at=yticks$y,labels=yticks$t)
    graphics::mtext(side=4,line=2,'t (Ma)')
    graphics::matlines(matrix(rep(xlim,nt),nrow=2),
                       rbind(lrmin,lrmax),lty=2,col='gray50')
}
