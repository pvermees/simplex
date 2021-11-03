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
#' st <- standard(preset='Temora',prefix=TEM)
#' dc <- drift(x=SHRIMP)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=st)
#' plot(cal,option=3)
#' }
#' @export
calibration <- function(lr,stand,snames=NULL,i=NULL,invert=FALSE,t=NULL){
    out <- lr
    out$standard <- stand
    dat <- subset(x=out,prefix=stand$prefix,snames=snames,i=i,invert=invert)
    if (stable(lr)) out$calibration <- stable_calibration(lr=dat)
    else out$calibration <- geochron_calibration(lr=dat,t=t)
    class(out) <- unique(append('calibration',class(out)))
    out
}

stable_calibration <- function(lr){
    LL <- function(par,lr){
        out <- 0
        np <- length(par)
        snames <- names(lr$samples)
        for (sname in snames){
            X <- lr$samples[[sname]]$lr$b0g[1:np]
            E <- lr$samples[[sname]]$lr$cov[1:np,1:np]
            LL <- stats::mahalanobis(x=X,center=par,cov=E)
            out <- out + LL
        }
        out
    }
    ni <- length(lr$method$num)
    init <- lr$samples[[1]]$lr$b0g[1:ni]
    wtdmean <- stats::optim(init,fn=LL,gr=NULL,method='BFGS',hessian=TRUE,lr=lr)
    out <- list()
    out$snames <- names(lr$samples)
    out$lr <- wtdmean$par
    out$cov <- solve(wtdmean$hessian)
    out
}

geochron_calibration <- function(lr,t=NULL,...){
    out <- list()
    oxide <- lr$method$oxide
    common <- do.call(lr$standard$fetchfun,args=list(dat=lr))$common
    type <- datatype(lr)
    if (type=="U-Pb"){
        out[[type]] <- geocal(lr,oxide=oxide['U'],t=t,type=type,
                              common=common['Pb206Pb204'])
    } else if (datatype(lr)=="U-Th-Pb"){
        out[['U-Pb']] <- geocal(lr,oxide=oxide['U'],t=t,
                                type="U-Pb",common=common['Pb206Pb204'])
        out[['Th-Pb']] <- geocal(lr,oxide=oxide['Th'],t=t,
                                 type="Th-Pb",common=common['Pb208Pb204'])
    } else {
        stop("Invalid data type.")
    }
    out
}

geocal <- function(lr,oxide,t=NULL,type,common){
    if (type=='U-Pb'){
        num <- c(oxide,'Pb206','Pb204')
        den <- c('U238','U238','Pb206')
    } else if (type=='Th-Pb'){
        num <- c(oxide,'Pb208','Pb204')
        den <- c('Th232','Th232','Pb208')
    } else {
        stop("Invalid data type.")
    }
    if (is.null(t)) t <- median(lr$samples[[1]]$time)
    out <- list(num=num,den=den,oxide=oxide,t=hours(t))
    out$common <- common
    out$snames <- names(lr$samples)
    yd <- beta2york(lr=lr,t=t,snames=out$snames,
                    num=num,den=den,common=common)
    out$fit <- IsoplotR::york(yd)
    out
}

beta2york <- function(lr,t=0,snames,num=c('UO2','Pb206','Pb204'),
                      den=c('U238','U238','Pb206'),common=0){
    if (missing(snames)) snames <- names(lr$samples)
    out <- matrix(0,nrow=length(snames),ncol=5)
    colnames(out) <- c('X','sX','Y','sY','rXY')
    rownames(out) <- snames
    for (sname in snames){
        LR <- lr$samples[[sname]]$lr
        b0g <- LR$b0g
        b0gnames <- names(b0g)
        nb0g <- length(b0gnames)
        OP <- paste0(num[1],'/',den[1])
        DP <- paste0(num[2],'/',den[2])
        ib0OP <- which(b0gnames %in% paste0('b0[',OP,']'))
        ib0DP <- which(b0gnames %in% paste0('b0[',DP,']'))
        igOP <- which(b0gnames %in% paste0('g[',OP,']'))
        igDP <- which(b0gnames %in% paste0('g[',DP,']'))
        bOP <- b0g[ib0OP] + b0g[igOP]*hours(t)
        bDP <- b0g[ib0DP] + b0g[igDP]*hours(t)
        X <- bOP
        Y <- bDP
        J <- matrix(0,2,length(LR$b0g))
        J[1,ib0OP] <- 1
        J[1,igOP] <- hours(t)
        J[2,ib0DP] <- 1
        J[2,igDP] <- hours(t)
        if (length(num)>2){ # for U-Pb and Th-Pb
            CD <- paste0(num[3],'/',den[3])
            ib0CD <- which(b0gnames %in% paste0('b0[',CD,']'))
            igCD <- which(b0gnames %in% paste0('g[',CD,']'))
            bCD <- b0g[ib0CD] + b0g[igCD]*hours(t)
            Y <- Y + log(1-exp(bCD)*common)
            J[2,ib0CD] <- -exp(bCD)*common/(1-exp(bCD)*common)
            J[2,igCD] <- -exp(bCD)*common*hours(t)/(1-exp(bCD)*common)
        }
        E <- J %*% LR$cov %*% t(J)
        out[sname,'X'] <- X
        out[sname,'Y'] <- Y
        out[sname,'sX'] <- sqrt(E[1,1])
        out[sname,'sY'] <- sqrt(E[2,2])
        out[sname,'rXY'] <- E[1,2]/(out[sname,'sX']*out[sname,'sY'])
    }
    out
}

#' @title plot calibration data
#' @description shows the calibration data on a logratio plot.
#' @param x an object of class \code{logratios}
#' @param type for \code{U-Th-Pb} data, if \code{type=2}, produces a
#'     Th-Pb calibration plot. Otherwise plots the U-Pb calibration.
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
plot.calibration <- function(x,option=1,...){
    if (stable(x)) {
        out <- calplot_stable(dat=x,...)
    } else {
        out <- calplot_geochronology(dat=x,option=option,...)
    }
    invisible(out)
}

calplot_stable <- function(dat,...){
    cal <- dat$calibration
    num <- dat$method$num
    den <- dat$method$den
    nn <- length(num)
    np <- nn*(nn-1)/2       # number of plot panels
    nc <- ceiling(sqrt(np)) # number of rows
    nr <- ceiling(np/nc)    # number of columns
    oldpar <- par(mfrow=c(nr,nc))
    ii <- 1
    snames <- cal$snames
    for (i in 1:(nn-1)){
        for (j in (i+1):nn){
            xlab <- paste0(num[i],'/',den[i])
            ylab <- paste0(num[j],'/',den[j])
            B <- beta2york(lr=dat,snames=snames,
                           num=num[c(i,j)],den=den[c(i,j)])
            IsoplotR::scatterplot(B,...)
            graphics::mtext(side=1,text=xlab,line=2)
            graphics::mtext(side=2,text=ylab,line=2)
            ell <- IsoplotR::ellipse(cal$lr[i],cal$lr[j],
                                     cal$cov[c(i,j),c(i,j)])
            graphics::polygon(ell,col='white')
            ii <- ii + 1
            if (ii>np) break
        }
    }
    par(oldpar)
}

calplot_geochronology <- function(dat,option=1,title=TRUE,...){
    cal <- dat$calibration
    np <- length(cal)       # number of plot panels
    nc <- ceiling(sqrt(np)) # number of columns
    nr <- ceiling(np/nc)    # number of rows
    oldpar <- par(mfrow=c(nr,nc))
    ii <- 1
    for (i1 in 1:(nr-1)){
        for (i2 in (i1+1):nc){
            if (ii>np) break
            calplot_geochronology_helper(dat,option=option,
                                         type=ii,title=title,...)
            ii <- ii + 1
        }
    }
}

calplot_geochronology_helper <- function(dat,option=1,type=1,title=TRUE,...){
    cal <- dat$calibration[[type]]
    num <- cal$num
    den <- cal$den
    yd <- beta2york(lr=dat,t=seconds(cal$t),snames=cal$snames,
                    num=num,den=den,common=cal$common)
    xlab <- paste0('log[',num[1],'/',den[1],']')
    ylab <- paste0('log[',num[2],'/',den[2],']')
    if (option==1){
        IsoplotR::scatterplot(yd,fit=cal$fit,xlab=xlab,ylab=ylab,...)
    } else {
        X <- NULL
        Y <- NULL
        for (sname in cal$snames){
            sp <- spot(dat=dat,sname=sname)
            Op <- betapars(spot=sp,ion=num[1])
            Pp <- betapars(spot=sp,ion=den[1])
            Dp <- betapars(spot=sp,ion=num[2])
            Cp <- betapars(spot=sp,ion=num[3])
            CD <- (Cp$sig-Cp$bkg)/(Dp$sig-Dp$bkg)
            CDdc <- exp(Dp$g*(Dp$t-Cp$t))
            newX <- log(Op$sig-Op$bkg) - log(Pp$sig-Pp$bkg) + Op$g*(Pp$t-Op$t)
            newY <- log(Dp$sig-Dp$bkg) - log(Pp$sig-Pp$bkg) + Dp$g*(Pp$t-Dp$t) +
                log(1 - cal$common*CD*CDdc)
            X <- cbind(X,newX)
            Y <- cbind(Y,newY)
        }
        xlim <- c(0,2)
        xlim[1] <- min(yd[cal$snames,'X']-3*yd[cal$snames,'sX'],X)
        xlim[2] <- max(yd[cal$snames,'X']+3*yd[cal$snames,'sX'],X)
        ylim <- c(0,2)
        ylim[1] <- min(yd[cal$snames,'Y']-3*yd[cal$snames,'sY'],Y)
        ylim[2] <- max(yd[cal$snames,'Y']+3*yd[cal$snames,'sY'],Y)
        IsoplotR::scatterplot(yd,fit=cal$fit,xlab=xlab,
                              ylab=ylab,xlim=xlim,ylim=ylim,...)
        graphics::matlines(X,Y,lty=1,col='darkgrey')
        if (option>2){
            graphics::points(X[1,],Y[1,],pch=21,bg='black')
            graphics::points(X[nrow(X),],Y[nrow(Y),],pch=21,bg='white')
        }
    }
    if (title) graphics::title(caltitle(fit=cal$fit))
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
