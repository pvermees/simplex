#' @title standard calibration
#' @description calibration of SIMS data using reference standards
#' @param lr an object of class \code{logratios}
#' @param stand an object of class \code{standard}
#' @param snames a vector with the sample names of selected standard
#'     analyses to be used for the calibration
#' @param i a vector of indices of selected standard analyses to be
#'     used for the calibration
#' @param invert if \code{TRUE}, inverts the selection made by
#'     \code{snames} or \code{i}
#' @param t analysis time that the signal should be regressed to
#' @return an object of class \code{calibration}
#' @examples
#' \dontrun{
#' data('SHRIMP',package='simplex')
#' st <- standard(preset='Temora',prefix='TEM')
#' dc <- drift(x=SHRIMP)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=st)
#' plot(cal,option=3)
#' }
#' @export
calibration <- function(lr,stand,pairing=NULL,prefix=NULL,
                        snames=NULL,i=NULL,invert=FALSE,t=NULL){
    out <- lr
    cal <- list()
    cal$stand <- stand
    cal$prefix <- prefix
    dat <- subset(x=out,prefix=prefix,snames=snames,i=i,invert=invert)
    if (is.null(t)) t <- stats::median(lr$samples[[1]]$time)
    cal$t <- t
    cal$tavg <- time_average(dat,t=t)
    if (is.null(pairing)){
        cal$cal <- average_calibration(tavg=cal$tavg,stand=stand)
    } else {
        cal$pairing <- pairing
        cal$cal <- regression_calibration(tavg=cal$tavg,
                                          pairing=pairing,stand=stand)
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
    LL <- function(par,tavg,i){
        out <- 0
        for (tav in tavg){
            out <- out +
                stats::mahalanobis(x=tav$val[i],center=par,cov=tav$cov[i,i])
        }
        out
    }
    ratios <- names(tavg[[1]]$val)
    i <- which(ratios %in% stand$ratios)
    wtdmean <- stats::optim(tavg[[1]]$val[i],fn=LL,gr=NULL,
                            method='BFGS',hessian=TRUE,tavg=tavg,i=i)
    val <- as.vector(wtdmean$par)
    covmat <- solve(wtdmean$hessian)
    list(ratios=ratios,val=unname(val),cov=unname(covmat))
}

regression_calibration <- function(tavg,pairing,stand){
    nr <- nrow(pairing)
    out <- as.data.frame(matrix(NA,nrow=nr,ncol=10))
    colnames(out) <- c('x','y','a','s[a]','b','s[b]','cov.ab','mswd','df','p.value')
    out[,1:2] <- pairing[,c('versus','smp')]
    for (i in 1:nr){
        yd <- beta2york_regression(tavg=tavg,pairing=pairing[i,],stand=stand)
        slope <- pairing[i,'slope']
        if (identical(slope,'auto')) fit <- IsoplotR::york(yd)
        else fit <- yorkfix(yd,b=as.numeric(slope))
        out[i,3:10] <- c(fit$a,fit$b,fit$cov.ab,fit$mswd,fit$df,fit$p.value)
    }
    out
}

beta2york_average <- function(tavg,i,j){
    snames <- names(tavg)
    out <- matrix(0,nrow=length(snames),ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        samp <- tavg[[sname]]
        out[sname,'X'] <- samp$val[i]
        out[sname,'Y'] <- samp$val[j]
        out[sname,'sX'] <- sqrt(samp$cov[i,i])
        out[sname,'sY'] <- sqrt(samp$cov[j,j])
        out[sname,'rXY'] <- samp$cov[i,j]/sqrt(samp$cov[i,i]*samp$cov[j,j])
    }
    out
}

beta2york_regression <- function(tavg,pairing,stand){
    snames <- names(tavg)
    out <- matrix(0,nrow=length(snames),ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    DP <- as.character(pairing[1,'smp'])
    OP <- as.character(pairing[1,'versus'])
    CD.smp <- as.character(pairing[1,'smp.c'])
    i.CD.smp <- match(pairing[1,'std.c'],stand$ratios)
    val.CD.std <- stand$val[i.CD.smp]
    J <- matrix(0,nrow=2,ncol=length(tavg[[1]]$val))
    colnames(J) <- names(tavg[[1]]$val)
    rownames(J) <- c('X','Y')
    J['X',OP] <- 1
    J['Y',DP] <- 1
    for (sname in snames){
        val <- tavg[[sname]]$val
        out[sname,'X'] <- val[OP]
        out[sname,'Y'] <- val[DP] - log(1 - exp(val[CD.smp]-val.CD.std))
        J['Y',CD.smp] <- -exp(val[CD.smp])/(1-exp(val[CD.smp]-val.CD.std))
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
#' @param x an object of class \code{logratios}
#' @param option if \code{option=1}, plots the best fit line through
#'     U-Pb and Th-Pb data. If \code{option=2}, adds the time-resolved
#'     raw data to the plot. If \code{option=3}, marks the first and
#'     last measurement by black and white circles, respectively.
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function.
#' @examples
#' \dontrun{
#' data('SHRIMP',package='simplex')
#' st <- standard(preset='Temora',prefix=TEM)
#' dc <- drift(x=SHRIMP)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=st)
#' plot(cal,option=3)
#' }
#' @method plot calibration
#'@export
plot.calibration <- function(dat,option=1,...){
    if (is.null(dat$calibration$pairing)){
        out <- calplot_stable(cal=dat$calibration,...)
    } else {
        out <- calplot_geochronology(dat=dat,option=option,...)
    }
    invisible(out)
}

calplot_stable <- function(cal,...){
    nn <- length(cal$cal$ratios)         # number of ratios
    snames <- names(cal$tavg)
    if (nn>1){
        np <- nn*(nn-1)/2       # number of panels
        nc <- ceiling(sqrt(np)) # number of rows
        nr <- ceiling(np/nc)    # number of columns
        oldpar <- graphics::par(mfrow=c(nr,nc))
        ratios <- cal$cal$ratios
        ij <- match(ratios,names(cal$tavg[[1]]$val))
        for (i in 1:(nn-1)){
            for (j in (i+1):nn){
                B <- beta2york_average(cal$tavg,ij[i],ij[j])
                IsoplotR::scatterplot(B,...)
                graphics::mtext(side=1,text=paste0('ln[',ratios[i],']'),line=2)
                graphics::mtext(side=2,text=paste0('ln[',ratios[j],']'),line=2)
                ell <- IsoplotR::ellipse(cal$cal$val[i],cal$cal$val[j],
                                         cal$cal$cov[c(i,j),c(i,j)])
                graphics::polygon(ell,col='white')
            }
        }
        graphics::par(oldpar)
    } else {
        ns <- length(snames)
        tab <- matrix(0,nrow=3,ncol=ns)
        rownames(tab) <- c('x','ll','ul')
        colnames(tab) <- snames
        tfact <- qnorm(0.975)
        for (sname in snames){
            tab['x',sname] <- cal$tavg[[sname]]$val
            dx <- tfact*sqrt(cal$tavg[[sname]]$cov)
            tab['ll',sname] <- tab['x',sname] - dx
            tab['ul',sname] <- tab['x',sname] + dx
        }
        matplot(rbind(1:ns,1:ns),tab[c('ll','ul'),],
                type='l',lty=1,col='black',bty='n',
                xlab='standard #',ylab='')
        points(1:ns,tab['x',],pch=16)
        lines(c(1,ns),rep(cal$cal$val,2),lty=2)
        lines(c(1,ns),rep(cal$cal$val-tfact*sqrt(cal$cal$cov),2),lty=3)
        lines(c(1,ns),rep(cal$cal$val+tfact*sqrt(cal$cal$cov),2),lty=3)
    }
}

calplot_geochronology <- function(dat=dat,option=option,...){
    cal <- dat$calibration
    nr <- nrow(cal$pairing)
    out <- as.data.frame(matrix(NA,nrow=nr,ncol=7))
    colnames(out) <- c('x','y','a','s[a]','b','s[b]','cov.ab')
    out[,1:2] <- cal$pairing[,c('versus','smp')]
    oldpar <- graphics::par(mfrow=c(1,nr))
    for (i in 1:nr){
        yd <- beta2york_regression(tavg=cal$tavg,
                                   pairing=cal$pairing[i,],
                                   stand=cal$stand)
        IsoplotR::scatterplot(yd,
                              xlab=paste0('ln[',cal$cal[i,'x'],']'),
                              ylab=paste0('ln[',cal$cal[i,'y'],']'))
        abline(a=cal$cal[i,'a'],b=cal$cal[i,'b'])
    }
    graphics::par(oldpar)
}

caltitle <- function(fit,sigdig=2,type=NA,...){
    args1 <- quote(a%+-%b~'(n='*n*')')
    args2 <- quote(a%+-%b)
    if (is.na(type)){
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
    }
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
