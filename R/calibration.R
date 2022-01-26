#' @title standard calibration
#' @description calibration of SIMS data using reference standards
#' @param lr an object of class \code{logratios}
#' @param stand an object of class \code{standard}
#' @param pairing (optional) an object of class \code{pairing}
#' @param prefix (optional) prefix of the aliquots to be used as
#'     calibration standards
#' @param snames (option) vector of sample names to be used as
#'     calibration standards
#' @param i (optional) vector of indices of aliquots to be used as
#'     calibration standards
#' @param invert if \code{TRUE}, inverts the selection made by
#'     \code{snames} or \code{i}
#' @param t analysis time that the signal should be regressed to
#' @return an object of class \code{calibration}
#' @examples
#' \dontrun{
#' data('SHRIMP_UPb',package='simplex')
#' st <- standard(preset='Temora')
#' dc <- drift(x=SHRIMP_UPb)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=st,prefix="TEM")
#' plot(cal)
#' }
#' @export
calibration <- function(lr,stand,pairing=NULL,prefix=NULL,
                        snames=NULL,i=NULL,invert=FALSE,t=NULL){
    out <- lr
    cal <- list()
    cal$stand <- stand
    cal$prefix <- prefix
    cal$snames <- subset2snames(dat=lr,prefix=prefix,snames=snames,i=i)
    dat <- subset(x=out,snames=cal$snames)
    if (is.null(t)) t <- stats::median(lr$samples[[1]]$time)
    cal$t <- t
    tavg <- time_average(dat,t=t)
    if (is.null(pairing)){
        cal$cal <- average_calibration(tavg=tavg,stand=stand)
    } else {
        cal$pairing <- pairing
        cal$cal <- regression_calibration(tavg=tavg,pairing=pairing,stand=stand)
    }
    out$calibration <- cal
    class(out) <- unique(append('calibration',class(out)))
    out
}

time_average <- function(lr,t=NULL){
    if (is.null(t)) t <- stats::median(lr$samples[[1]]$time)
    tt <- hours(t)
    snames <- names(lr$samples)
    rnames <- paste0(lr$method$num,'/',lr$method$den)
    out <- list()
    nlr <- length(lr$samples[[1]]$lr$b0g)/2
    J <- cbind(diag(nlr),tt*diag(nlr))
    rownames(J) <- rnames
    ib0 <- 1:nlr
    ig <- (nlr+1):(2*nlr)
    for (sname in snames){
        LR <- lr$samples[[sname]]$lr
        out[[sname]] <- list()
        out[[sname]]$val <- as.vector(c(1,tt) %*% rbind(LR$b0g[ib0],LR$b0g[ig]))
        names(out[[sname]]$val) <- rnames
        out[[sname]]$cov <- J %*% LR$cov %*% t(J)
    }
    out
}

average_calibration <- function(tavg,stand){
    LL <- function(par,tavg,ratios){
        out <- 0
        for (tav in tavg){
            out <- out +
                stats::mahalanobis(x=tav$val[ratios],center=par,
                                   cov=tav$cov[ratios,ratios])
        }
        out
    }
    ratios <- intersect(names(tavg[[1]]$val),names(stand$val))
    wtdmean <- stats::optim(tavg[[1]]$val[ratios],fn=LL,gr=NULL,
                            method='BFGS',hessian=TRUE,
                            tavg=tavg,ratios=ratios)
    list(val=wtdmean$par,cov=solve(wtdmean$hessian))
}

regression_calibration <- function(tavg,pairing,stand){
    nr <- nrow(pairing)
    cnames <- c('a','s[a]','b','s[b]','cov.ab','mswd','df','p.value')
    nc <- length(cnames)
    out <- as.data.frame(matrix(NA,nrow=nr,ncol=nc))
    colnames(out) <- cnames
    for (i in 1:nr){
        yd <- beta2york_regression(tavg=tavg,pairing=pairing[i,],stand=stand)
        slope <- pairing[i,'slope']
        if (identical(slope,'auto')) fit <- IsoplotR::york(yd)
        else fit <- yorkfix(yd,b=as.numeric(slope))
        out[i,] <- c(fit$a,fit$b,fit$cov.ab,fit$mswd,fit$df,fit$p.value)
    }
    out
}

beta2york_average <- function(tavg,xratio,yratio){
    snames <- names(tavg)
    out <- matrix(0,nrow=length(snames),ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        samp <- tavg[[sname]]
        out[sname,'X'] <- samp$val[xratio]
        out[sname,'Y'] <- samp$val[yratio]
        vX <- samp$cov[xratio,xratio]
        vY <- samp$cov[yratio,yratio]
        sXY <- samp$cov[xratio,yratio]
        out[sname,'sX'] <- sqrt(vX)
        out[sname,'sY'] <- sqrt(vY)
        out[sname,'rXY'] <- sXY/sqrt(vX*vY)
    }
    out
}

beta2york_regression <- function(tavg,pairing,stand){
    snames <- names(tavg)
    out <- matrix(0,nrow=length(snames),ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    J <- matrix(0,nrow=2,ncol=length(tavg[[1]]$val))
    colnames(J) <- names(tavg[[1]]$val)
    rownames(J) <- c('X','Y')
    J['X',pairing$X] <- 1
    J['Y',pairing$Y] <- 1
    for (sname in snames){
        val <- tavg[[sname]]$val
        if (stand$measured){
            CD <- 0
        } else {
            CD <- log(1-exp(val[pairing$C]-stand$val[pairing$C]))
            J['Y',pairing$C] <-
                -exp(val[pairing$C])/(1-exp(val[pairing$C]-stand$val[pairing$C]))
        }
        out[sname,'X'] <- val[pairing$X]
        out[sname,'Y'] <- val[pairing$Y] - CD
        E <- J %*% tavg[[sname]]$cov %*% t(J)
        out[sname,'sX'] <- sqrt(E['X','X'])
        out[sname,'sY'] <- sqrt(E['Y','Y'])
        out[sname,'rXY'] <- E['X','Y']/sqrt(E['X','X']*E['Y','Y'])
    }
    out
}

yorkfix <- function(xy,b,alpha=0.05){
    SS <- function(a,b,xy){
        XY <- IsoplotR:::get.york.xy(xy,a,b)
        dX <- XY[,1]-xy[,'X']
        sX <- xy[,'sX']
        dY <- XY[,2]-xy[,'Y']
        sY <- xy[,'sY']
        rXY <- xy[,'rXY']
        SS <- (dX/sX)^2 + (dY/sY)^2 - 2*rXY*dX*dY/(sX*sY)
        sum(SS)
    }
    init <- lm(Y - b*X ~ 1, data=as.data.frame(xy))$coefficients
    fit <- stats::optim(init,SS,method='BFGS',b=b,xy=xy,hessian=TRUE)
    out <- list(type='york')
    a <- fit$par
    sa <- sqrt(solve(fit$hessian))
    out$a <- c(a,sa)
    names(out$a) <- c('a','s[a]')
    out$b <- c(b,0)
    names(out$b) <- c('b','s[b]')
    out$cov.ab <- 0
    out$df <- nrow(xy)-1
    out$mswd <- fit$value/out$df
    out$p.value <- as.numeric(1-stats::pchisq(fit$value,out$df))
    out$alpha <- alpha
    out
}

#' @title plot calibration data
#' @description shows the calibration data on a logratio plot.
#' @param dat an object of class \code{logratios}
#' @param show.numbers logical. If \code{TRUE}, labels the error
#'     ellipses with the aliquot numbers.
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function.
#' @examples
#' \dontrun{
#' data('SHRIMP_UPb',package='simplex')
#' st <- standard(preset='Temora')
#' dc <- drift(x=SHRIMP_UPb)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=st,prefix='TEM')
#' plot(cal)
#' }
#' @method plot calibration
#'@export
plot.calibration <- function(dat,show.numbers=TRUE,...){
    if (is.null(dat$calibration$pairing)){
        out <- calplot_stable(dat=dat,show.numbers=show.numbers,...)
    } else {
        out <- calplot_geochronology(dat=dat,show.numbers=show.numbers,...)
    }
    invisible(out)
}

calplot_stable <- function(dat,show.numbers=TRUE,...){
    cal <- dat$calibration$cal
    ratios <- names(cal$val)
    nrat <- length(ratios)
    caldat <- subset(x=dat,snames=dat$calibration$snames)
    tavg <- time_average(caldat,t=dat$calibration$t)
    if (nrat>1){
        np <- nrat*(nrat-1)/2   # number of panels
        nc <- ceiling(sqrt(np)) # number of rows
        nr <- ceiling(np/nc)    # number of columns
        oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,1,1))
        for (i in 1:(nrat-1)){
            for (j in (i+1):nrat){
                B <- beta2york_average(tavg,ratios[i],ratios[j])
                IsoplotR::scatterplot(B,...)
                if (show.numbers){
                    istd <- which(names(dat$samples) %in% dat$calibration$snames)
                    text(x=B[,'X'],y=B[,'Y'],labels=istd,pos=3,offset=0.1)
                }
                graphics::mtext(side=1,text=paste0('ln[',ratios[i],']'),line=2)
                graphics::mtext(side=2,text=paste0('ln[',ratios[j],']'),line=2)
                ell <- IsoplotR::ellipse(cal$val[i],cal$val[j],
                                         cal$cov[c(i,j),c(i,j)])
                graphics::polygon(ell,col='white')
            }
        }
        graphics::par(oldpar)
    } else {
        snames <- dat$calibration$snames
        ns <- length(snames)
        tab <- matrix(0,nrow=3,ncol=ns)
        rownames(tab) <- c('x','ll','ul')
        colnames(tab) <- snames
        tfact <- qnorm(0.975)
        for (sname in snames){
            tab['x',sname] <- tavg[[sname]]$val
            dx <- tfact*sqrt(tavg[[sname]]$cov)
            tab['ll',sname] <- tab['x',sname] - dx
            tab['ul',sname] <- tab['x',sname] + dx
        }
        matplot(rbind(1:ns,1:ns),tab[c('ll','ul'),],
                type='l',lty=1,col='black',bty='n',
                xlab='standard #',ylab='')
        dsd <- IsoplotR:::roundit(cal$val,tfact*sqrt(cal$cov),sigdig=2)
        del <- dsd[1]
        err <- dsd[2]
        rat <- names(cal$val)
        tit <- substitute(paste(rat,'=',del %+-% err,' (95% CI)'))
        mtext(text=tit)
        points(1:ns,tab['x',],pch=16)
        lines(c(1,ns),rep(cal$val,2),lty=2)
        lines(c(1,ns),rep(cal$val-tfact*sqrt(cal$cov),2),lty=3)
        lines(c(1,ns),rep(cal$val+tfact*sqrt(cal$cov),2),lty=3)
    }
}

calplot_geochronology <- function(dat=dat,option=option,show.numbers=TRUE,...){
    cal <- dat$calibration$cal
    pair <- dat$calibration$pairing
    nr <- nrow(pair)
    cnames <- c('a','s[a]','b','s[b]','cov.ab')
    nc <- length(cnames)
    out <- as.data.frame(matrix(NA,nrow=nr,ncol=nc))
    oldpar <- graphics::par(mfrow=c(1,nr),mar=c(3.5,3.5,3.5,1),mgp=c(2,0.5,0))
    snames <- dat$calibration$snames
    tavg <- time_average(subset(x=dat,snames=snames),t=dat$calibration$t)
    for (i in 1:nr){
        yd <- beta2york_regression(tavg=tavg,
                                   pairing=pair[i,],
                                   stand=dat$calibration$stand)
        fit <- cal2york(cal[i,])
        IsoplotR::scatterplot(yd,fit=fit,
                              xlab=paste0('ln[',pair[i,'X'],']'),
                              ylab=paste0('ln[',pair[i,'Y'],']'))
        if (show.numbers){
            istd <- which(names(dat$samples) %in% snames)
            text(x=yd[,'X'],y=yd[,'Y'],labels=istd,pos=3,offset=0.1)
        }
        caltitle(fit,...)
    }
    graphics::par(oldpar)
}

cal2york <- function(cal){
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

caltitle <- function(fit,sigdig=2,...){
    args1 <- quote(a%+-%b~'(n='*n*')')
    args2 <- quote(a%+-%b)
    intercept <- IsoplotR:::roundit(fit$a[1],fit$a[2],sigdig=sigdig)
    slope <- IsoplotR:::roundit(fit$b[1],fit$b[2],sigdig=sigdig)
    expr1 <- 'slope ='
    expr2 <- 'intercept ='
    list1 <- list(a=slope[1],
                  b=slope[2],
                  u='',
                  n=fit$df+2)
    list2 <- list(a=intercept[1],
                  b=intercept[2],
                  u='')
    call1 <- substitute(e~a,list(e=expr1,a=args1))
    call2 <- substitute(e~a,list(e=expr2,a=args2))
    line1 <- do.call(substitute,list(eval(call1),list1))
    line2 <- do.call(substitute,list(eval(call2),list2))
    line3 <- substitute('MSWD ='~a*', p('*chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    IsoplotR:::mymtext(line1,line=2,...)
    IsoplotR:::mymtext(line2,line=1,...)
    IsoplotR:::mymtext(line3,line=0,...)
}
