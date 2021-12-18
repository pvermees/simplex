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
    if (is.null(cal$calibration$pairing)){
        out <- calibrate_average(dat=cal,exterr=exterr)
    } else {
        out <- calibrate_regression(dat=cal,exterr=exterr)
    }
    out
}

calibrate_average <- function(dat,exterr=FALSE){
    out <- dat
    cal <- dat$calibration$cal
    stand <- dat$calibration$stand
    ns <- length(dat$samples)
    tavg <- time_average(dat,t=dat$calibration$t)
    alliso <- names(tavg[[1]]$val)
    staiso <- names(stand$val)
    caliso <- names(cal$val) # caliso is a subset of staiso
    nai <- length(alliso)
    nsi <- length(staiso)
    nci <- length(caliso)
    E <- matrix(0,nrow=ns*nai+nsi+nci,ncol=ns*nai+nsi+nci)
    J <- matrix(0,nrow=ns*nai,ncol=ns*nai+nsi+nci)
    ical <- which(alliso %in% caliso)
    ista <- which(staiso %in% caliso)
    E[ns*nai+1:nsi,ns*nai+1:nsi] <- as.matrix(stand$cov)
    E[ns*nai+nsi+1:nci,ns*nai+nsi+1:nci] <- as.matrix(cal$cov)
    val <- rep(0,ns*nai)
    for (i in 1:ns){
        iall <- (i-1)*nai+1:nai
        val[iall] <- tavg[[i]]$val
        E[iall,iall] <- tavg[[i]]$cov
        J[iall,iall] <- diag(nai)
        val[iall[ical]] <- val[iall[ical]] + stand$val[ista] - cal$val
        if (exterr){
            J[iall[ista],ns*nai+1:nsi] <- diag(nsi)[ista,]
            J[iall[ical],ns*nai+nsi+1:nci] <- -diag(nci)
        }
    }
    out$calibrated <- list(snames=names(dat$samples),ratios=alliso,
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
    E[istd,istd] <- as.matrix(cal$stand$cov)
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
                cal$stand$val[iDPstd[j]]
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
plot.calibrated <- function(x,...){
    if (is.null(x$calibration$pairing)){
        out <- caldplot_stable(dat=x,...)
    } else {
        out <- caldplot_geochronology(dat=x,option=option,...)
    }
    invisible(out)
}

caldplot_stable <- function(dat,...){
    sta <- dat$calibration$stand
    cal <- dat$calibrated
    tab <- data2table.calibrated(dat)
    nrat <- length(cal$ratios)
    if (nrat>1){
        np <- nrat*(nrat-1)/2   # number of plot panels
        nc <- ceiling(sqrt(np)) # number of rows
        nr <- ceiling(np/nc)    # number of columns
        oldpar <- graphics::par(mfrow=c(nr,nc),mar=rep(3.5,4))
        for (i in 1:(nrat-1)){
            for (j in (i+1):nrat){
                xratio <- cal$ratios[i]
                yratio <- cal$ratios[j]
                Xc <- tab[,xratio]
                Yc <- tab[,yratio]
                sXc <- tab[,paste0('s[',xratio,']')]
                sYc <- tab[,paste0('s[',yratio,']')]
                rXYc <- tab[,paste0('r[',xratio,',',yratio,']')]
                Xs <- sta$val[xratio]
                sXs <- sqrt(sta$cov[xratio,xratio])
                Ys <- sta$val[yratio]
                sYs <- sqrt(sta$cov[yratio,yratio])
                xlim <- c(min(Xc-3*sXc,Xs-3*sXs),max(Xc+3*sXc,Xs+3*sXs))
                ylim <- c(min(Yc-3*sYc,Ys-3*sYs),max(Yc+3*sYc,Ys+3*sYs))
                graphics::plot(xlim,ylim,type='n',ann=FALSE)
                deltagrid(dat,xratio,yratio)
                IsoplotR::scatterplot(cbind(Xc,sXc,Yc,sYc,rXYc),
                                      xlim=xlim,ylim=ylim,add=TRUE,...)
                graphics::mtext(side=1,text=paste0('ln[',xratio,']'),line=2)
                graphics::mtext(side=2,text=paste0('ln[',yratio,']'),line=2)
            }
        }
    } else {
        oldpar <- graphics::par(mar=c(3.5,3.5,1.5,3.5),mgp=c(2,0.75,0))
        ylab <- cal$ratios
        ns <- length(cal$snames)
        tfact <- qnorm(0.975)
        lr <- tab[,1]
        ll <- lr - tfact*tab[,2]
        ul <- lr + tfact*tab[,2]
        Ys <- sta$val[ylab]
        sYs <- as.numeric(sqrt(sta$cov[ylab,ylab]))
        xlim <- c(1,ns)
        ylim <- c(min(Ys-3*sYs,ll),max(Ys+3*sYs,ul))
        graphics::plot(xlim,ylim,type='n',xlab='sample #',ylab=ylab)
        deltagrid(dat)
        matlines(rbind(1:ns,1:ns),rbind(ll,ul),lty=1,col='black')
        points(1:ns,lr,pch=16)
        xlim <- graphics::par('usr')[1:2]
        lines(xlim,rep(Ys,2),lty=1,col='red')
        lines(xlim,rep(Ys-tfact*sYs,2),lty=2)
        lines(xlim,rep(Ys+tfact*sYs,2),lty=2)
    }
    graphics::par(oldpar)
}

deltagrid <- function(dat,xratio,yratio){
    if (missing(xratio) | missing(yratio)){
        multipanel <- FALSE
        xratio <- dat$calibrated$ratios
        yratio <- xratio
    } else {
        multipanel <- TRUE
    }
    usr <- graphics::par('usr')
    xlim <- usr[1:2]
    ylim <- usr[3:4]
    sta <- dat$calibration$stand
    staratios <- names(sta$val)
    calratios <- names(dat$calibrated$val)
    if (multipanel){
        x <- sta$val[xratio]
        dxmin <- 1000*(xlim[1]-sta$ref$val[xratio])
        dxmax <- 1000*(xlim[2]-sta$ref$val[xratio])
        dxticks <- pretty(c(dxmin,dxmax))
        xticks <- sta$ref$val[xratio] + dxticks/1000
        nxt <- length(xticks)
    }
    y <- sta$val[yratio]
    dymin <- 1000*(ylim[1]-sta$ref$val[yratio])
    dymax <- 1000*(ylim[2]-sta$ref$val[yratio])
    dyticks <- pretty(c(dymin,dymax))
    yticks <- sta$ref$val[yratio] + dyticks/1000
    nyt <- length(yticks)
    if (multipanel){
        graphics::matlines(rbind(xticks,xticks),
                           matrix(rep(ylim,nxt),ncol=nxt),
                           lty=3,col='black')
        graphics::matlines(matrix(rep(xlim,nyt),ncol=nyt),
                           rbind(yticks,yticks),lty=3,col='black')
        graphics::axis(side=3,at=xticks,labels=dxticks)
        graphics::mtext(expression(delta*"'"),side=3,line=2)
        graphics::lines(rep(x,2),ylim,lty=2,col='red')
        graphics::text(x,ylim[1],labels=dat$calibration$prefix,
                       pos=4,srt=90,offset=0)
        graphics::lines(xlim,rep(y,2),lty=2,col='red')
        graphics::text(xlim[1],y,pos=4,offset=0,
                       labels=dat$calibration$prefix)
    } else {
        graphics::matlines(matrix(rep(xlim,nyt),ncol=nyt),
                           rbind(yticks,yticks),lty=3,col='black')
    }
    graphics::axis(side=4,at=yticks,labels=dyticks)
    graphics::mtext(expression(delta*"'"),side=4,line=2)
}

caldplot_geochronology <- function(dat,option=1,...){
    cal <- dat$calibration$cal
    pairing <- dat$calibration$pairing
    stand <- dat$calibration$stand
    tab <- data2table.logratios(dat,t=dat$calibration$t)
    nr <- nrow(pairing)
    oldpar <- graphics::par(mfrow=c(1,nr))
    for (i in 1:nr){
        X <- paste0('ln[',pairing[i,'versus'],']')
        Y <- paste0('ln[',pairing[i,'smp'],']')
        sX <- paste0('s(',X,')')
        sY <- paste0('s(',Y,')')
        rXY <- paste0('r(',X,',',Y,')')
        if (!(rXY %in% colnames(tab)))
            rXY <- paste0('r(',Y,',',X,')')
        XY <- tab[,c(X,sX,Y,sY,rXY)]
        xlim <- c(min(XY[,1]-3*XY[,2]),max(XY[,1]+3*XY[,2]))
        ylim <- c(min(XY[,3]-3*XY[,4]),max(XY[,3]+3*XY[,4]))
        fit <- cald2york(cal[i,])
        plot(xlim,ylim,type='n',xlab=X,ylab=Y)
        agegrid(fit=fit,pairing=pairing[i,],stand=stand)
        IsoplotR::scatterplot(XY,fit=fit,add=TRUE)
    }
    graphics::par(oldpar)
}

agegrid <- function(fit,pairing,stand){
    if (pairing$std=='Pb206/U238'){
        lambda <- IsoplotR::settings('lambda','U238')[1]
    } else if (pairing$std=='Pb208/Th232'){
        lambda <- IsoplotR::settings('lambda','Th232')[1]
    } else {
        return(NA)
    }
    usr <- graphics::par('usr')
    xlim <- usr[1:2]
    ylim <- usr[3:4]
    a <- fit$a[1]
    b <- fit$b[1]
    yrange <- c(ylim[1] - b*diff(xlim),
                ylim[2] + b*diff(xlim))
    logDPstd <- stand$val[match(pairing$std,stand$ratios)]
    logDPlim <- logDPstd + yrange - a - b * xlim
    tlim <- log(exp(logDPlim)+1)/lambda
    tticks <- pretty(tlim)
    logDPticks <- log(exp(lambda*tticks)-1)
    nt <- length(tticks)
    xl <- rep(xlim[1],nt)
    xu <- rep(xlim[2],nt)
    yl <- logDPticks - logDPstd + a + b * xl
    yu <- logDPticks - logDPstd + a + b * xu
    matlines(rbind(xl,xu),rbind(yl,yu),lty=3,col='black')
    top <- (yu>ylim[2])
    axis(side=3,at=xlim[1]+(ylim[2]-yl[top])/b,labels=tticks[top])
    axis(side=4,at=yu[!top],labels=tticks[!top])
}

cald2york <- function(cal){
    out <- list()
    out$a <- as.numeric(cal[c('a','s[a]')])
    out$b <- as.numeric(cal[c('b','s[b]')])
    out$cov.ab <- as.numeric(cal['cov.ab'])
    out$df <- as.numeric(cal['df'])
    out$mswd <- as.numeric(cal['mswd'])
    out$p.value <- as.numeric(cal['p.value'])
    out$type <- 'york'
    out$alpha <- 0.05
    out
}
