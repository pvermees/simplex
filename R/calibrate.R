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
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28-O'))
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
    std <- dat$calibration$stand
    ns <- length(dat$samples)
    tavg <- time_average(dat,t=dat$calibration$t)
    alliso <- names(tavg[[1]]$val)
    stdiso <- names(std$val)
    caliso <- names(cal$val) # caliso is the intersection of alliso and stdiso
    nai <- length(alliso)
    nsi <- length(stdiso)
    nci <- length(caliso)
    E <- matrix(0,nrow=ns*nai+nsi+nci,ncol=ns*nai+nsi+nci)
    J <- matrix(0,nrow=ns*nai,ncol=ns*nai+nsi+nci)
    iallincal <- match(caliso,alliso)
    istdincal <- match(caliso,stdiso)
    E[ns*nai+1:nsi,ns*nai+1:nsi] <- as.matrix(std$cov)
    E[ns*nai+nsi+1:nci,ns*nai+nsi+1:nci] <- as.matrix(cal$cov)
    val <- rep(0,ns*nai)
    for (i in 1:ns){
        iall <- (i-1)*nai+1:nai
        val[iall] <- tavg[[i]]$val
        E[iall,iall] <- tavg[[i]]$cov
        J[iall,iall] <- diag(nai)
        val[iall[iallincal]] <- val[iall[iallincal]] + std$val[istdincal] - cal$val
        if (exterr){
            J[iall[iallincal],ns*nai+1:nsi] <- diag(nsi)[istdincal,]
            J[iall[iallincal],ns*nai+nsi+1:nci] <- -diag(nci)
        }
    }
    out$calibrated <- list(snames=names(dat$samples),ratios=alliso,
                           val=val,cov=J %*% E %*% t(J))
    class(out) <- unique(append("calibrated",class(out)))
    out
}

calibrate_regression <- function(dat,exterr=FALSE){
    out <- dat
    cal <- dat$calibration$cal
    std <- dat$calibration$stand
    pairing <- dat$calibration$pairing
    tavg <- time_average(dat,t=dat$calibration$t)
    alliso <- names(tavg[[1]]$val)
    outiso <- alliso[!(alliso %in% pairing$X)]
    stdiso <- names(std$val)
    ns <- length(tavg)
    nai <- length(alliso)
    noi <- length(outiso)
    nsi <- length(stdiso)
    nab <- nrow(pairing)
    iout <- which(alliso %in% outiso)
    istd <- ns*nai+1:nsi
    iab <- matrix(0,nab,2)
    iDPa <- rep(NA,nab)
    iDPo <- rep(NA,nab)
    iDPs <- rep(NA,nab)
    iOP <- rep(NA,nab)
    for (i in 1:nab){
        iab[i,] <- ns*nai+nsi+(i-1)*2+1:2
        iDPa[i] <- which(alliso %in% pairing[i,'Y'])
        iDPo[i] <- which(outiso %in% pairing[i,'Y'])
        iDPs[i] <- which(stdiso %in% pairing[i,'Y'])
        iOP[i] <- which(alliso %in% pairing[i,'X'])
    }
    val <- rep(0,ns*noi)
    E <- matrix(0,nrow=ns*nai+nsi+2*nab,ncol=ns*nai+nsi+2*nab)
    E[istd,istd] <- as.matrix(std$cov)
    for (i in 1:nab){
        E[iab[i,1],iab[i,1]] <- cal[i,'s[a]']^2
        E[iab[i,2],iab[i,2]] <- cal[i,'s[b]']^2
        E[iab[i,1],iab[i,2]] <- cal[i,'cov.ab']
        E[iab[i,2],iab[i,1]] <- cal[i,'cov.ab']
    }
    J <- matrix(0,nrow=ns*noi,ncol=ns*nai+nsi+2*nab)
    for (i in 1:ns){
        ioi <- (i-1)*noi+1:noi
        val[ioi] <- tavg[[i]]$val[iout]
        iai <- (i-1)*nai+1:nai
        E[iai,iai] <- tavg[[i]]$cov
        J[ioi,iai] <- diag(nai)[iout,] # dval/dval
        for (j in 1:nab){
            val[ioi[iDPo[j]]] <- tavg[[i]]$val[iDPa[j]] -
                cal[j,'a'] - cal[j,'b']*tavg[[i]]$val[iOP[j]] +
                std$val[iDPs[j]]
            J[ioi[iDPo[j]],iai[iOP[j]]] <- -cal[j,'b'] # dval/dOP
            if (exterr){
                J[ioi[iDPo[j]],iab[j,1]] <- -1 # dval/da
                J[ioi[iDPo[j]],iab[j,2]] <- -tavg[[i]]$val[iOP[j]] # dval/db
                J[ioi[iDPo[j]],istd[iDPs[j]]] <- 1 # dval/dDPs
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
#' @param show.numbers logical. If \code{TRUE}, numbers the error
#'     ellipses by aliquot.
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function
#' @examples
#' \dontrun{
#' data('Cameca_UPb',package='simplex')
#' dc <- drift(x=Cameca_UPb)
#' lr <- logratios(x=dc)
#' st <- standard(preset='Temora-t')
#' p <- pairing(lr,stand=st)
#' cal <- calibration(lr=lr,pairing=p,stand=st,prefix='Tem')
#' cd <- calibrate(cal)
#' plot(cd)
#' }
#' @method plot calibrated
#' @export
plot.calibrated <- function(x,show.numbers=TRUE,...){
    if (is.null(x$calibration$pairing)){
        out <- caldplot_stable(dat=x,show.numbers=show.numbers,...)
    } else {
        out <- caldplot_geochronology(dat=x,show.numbers=show.numbers,...)
    }
    invisible(out)
}

caldplot_stable <- function(dat,calfit=FALSE,show.numbers=TRUE,...){
    sta <- dat$calibration$stand
    cal <- dat$calibrated
    tab <- data2table.calibrated(dat,log4lab=FALSE)
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
                xlim <- c(min(Xc-3*sXc),max(Xc+3*sXc))
                ylim <- c(min(Yc-3*sYc),max(Yc+3*sYc))
                if (calfit){
                    Xs <- sta$val[xratio]
                    sXs <- sqrt(sta$cov[xratio,xratio])
                    Ys <- sta$val[yratio]
                    sYs <- sqrt(sta$cov[yratio,yratio])
                    xlim[1] <- min(xlim[1],min(Xs-3*sXs))
                    xlim[2] <- max(xlim[2],max(Xs+3*sXs))
                    ylim[1] <- min(ylim[1],min(Ys-3*sYs))
                    ylim[2] <- max(ylim[2],max(Ys+3*sYs))
                }
                graphics::plot(xlim,ylim,type='n',ann=FALSE)
                if (!is.null(sta$ref)) deltagrid(dat,xratio,yratio)
                IsoplotR::scatterplot(cbind(Xc,sXc,Yc,sYc,rXYc),
                                      xlim=xlim,ylim=ylim,add=TRUE,
                                      show.numbers=show.numbers,...)
                graphics::mtext(side=1,text=paste0('ln[',xratio,']'),line=2)
                graphics::mtext(side=2,text=paste0('ln[',yratio,']'),line=2)
            }
        }
    } else {
        oldpar <- graphics::par(mar=c(3.5,3.5,1.5,3.5),mgp=c(2,0.75,0))
        ylab <- cal$ratios
        ns <- length(cal$snames)
        tfact <- stats::qnorm(0.975)
        lr <- tab[,1]
        ll <- lr - tfact*tab[,2]
        ul <- lr + tfact*tab[,2]
        Ys <- sta$val[ylab]
        sYs <- as.numeric(sqrt(sta$cov[ylab,ylab]))
        xlim <- c(1,ns)
        ylim <- c(min(Ys-3*sYs,ll),max(Ys+3*sYs,ul))
        graphics::plot(xlim,ylim,type='n',xlab='sample #',ylab=ylab)
        if (!null(sta$ref)) deltagrid(dat)
        graphics::matlines(rbind(1:ns,1:ns),rbind(ll,ul),lty=1,col='black')
        graphics::points(1:ns,lr,pch=16)
        xlim <- graphics::par('usr')[1:2]
        graphics::lines(xlim,rep(Ys,2),lty=1,col='red')
        graphics::lines(xlim,rep(Ys-tfact*sYs,2),lty=2)
        graphics::lines(xlim,rep(Ys+tfact*sYs,2),lty=2)
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

caldplot_geochronology <- function(dat,calfit=FALSE,show.numbers=TRUE,...){
    cal <- dat$calibration$cal
    pairing <- dat$calibration$pairing
    stand <- dat$calibration$stand
    tab <- data2table.logratios(dat,t=dat$calibration$t)
    nr <- nrow(pairing)
    oldpar <- graphics::par(mfrow=c(1,nr),mar=rep(3,4),mgp=c(1.5,0.5,0))
    for (i in 1:nr){
        X <- paste0('ln[',pairing[i,'X'],']')
        Y <- paste0('ln[',pairing[i,'Y'],']')
        sX <- paste0('s(',X,')')
        sY <- paste0('s(',Y,')')
        rXY <- paste0('r(',X,',',Y,')')
        if (!(rXY %in% colnames(tab))) rXY <- paste0('r(',Y,',',X,')')
        xy <- tab[,c(X,sX,Y,sY,rXY)]
        xlim <- c(min(xy[,1]-3*xy[,2]),max(xy[,1]+3*xy[,2]))
        ylim <- c(min(xy[,3]-3*xy[,4]),max(xy[,3]+3*xy[,4]))
        if (calfit){
            yc <- cal[i,'a'] + cal[i,'b']*tab[,X]
            ylim[1] <- min(ylim[1],min(yc))
            ylim[2] <- max(ylim[2],max(yc))
        }
        fit <- cal2york(cal[i,])
        graphics::plot(xlim,ylim,type='n',xlab=X,ylab=Y,...)
        agegrid(fit=fit,pairing=pairing[i,],stand=stand)
        IsoplotR::scatterplot(xy,fit=fit,add=TRUE,show.numbers=show.numbers)
    }
    graphics::par(oldpar)
}

agegrid <- function(fit,pairing,stand){
    if (identical(pairing$Y,'Pb206/U238')){
        lambda <- IsoplotR::settings('lambda','U238')[1]
    } else if (identical(pairing$Y,'Pb208/Th232')){
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
    logDPstd <- stand$val[pairing$Y]
    logDPlim <- logDPstd + yrange - a - b * xlim
    tlim <- log(exp(logDPlim)+1)/lambda
    tticks <- pretty(tlim)
    logDPticks <- log(exp(lambda*tticks)-1)
    nt <- length(tticks)
    xl <- rep(xlim[1],nt)
    xu <- rep(xlim[2],nt)
    yl <- logDPticks - logDPstd + a + b * xl
    yu <- logDPticks - logDPstd + a + b * xu
    graphics::matlines(rbind(xl,xu),rbind(yl,yu),lty=3,col='black')
    top <- (yu>ylim[2])
    graphics::axis(side=3,at=xlim[1]+(ylim[2]-yl[top])/b,labels=tticks[top])
    graphics::axis(side=4,at=yu[!top],labels=tticks[!top])
    graphics::mtext('age (Ma)',side=3,line=1.5)
    graphics::mtext('age (Ma)',side=4,line=1.5)
}
