#' @title calibrate SIMS data
#' @description convert signal logratios to atomic logratios
#' @param cal an object of class \code{calibration}
#' @param exterr logical flag indicating whether the uncertainty of
#'     the calibration should be propagated.
#' @return an object of class \code{calibrated}
#' @examples
#' data('oxygen')
#' dc <- drift(x=oxygen)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='NBS28'))
#' cd <- calibrate(cal)
#' del <- delta(cd,log=FALSE)
#' tab <- data2table(del)
#' @export
calibrate <- function(cal,exterr=FALSE){
    if (stable(cal)) out <- calibrate_stable(dat=cal,exterr=exterr)
    else out <- calibrate_geochron(dat=cal,exterr=exterr)
    out
}

calibrate_stable <- function(dat,exterr=FALSE){
    out <- dat
    scal <- do.call(dat$standard$fetchfun,args=list(dat=dat))
    dcal <- dat$calibration
    num <- dat$m$num
    den <- dat$m$den
    snames <- names(dat$samples)
    ns <- length(snames)
    nr <- length(num)
    calibrated <- list()
    calibrated$num <- num
    calibrated$den <- den
    calibrated$lr <- rep(0,nr*ns)
    E <- matrix(0,nr*(ns+2),nr*(ns+2))
    J <- matrix(0,nr*ns,nr*(ns+2))
    J[1:(ns*nr),1:(ns*nr)] <- diag(ns*nr)
    for (i in 1:ns){
        sname <- snames[i]
        selected <- (i-1)*nr + (1:nr)
        sp <- spot(dat,sname=sname)
        calibrated$lr[selected] <- sp$lr$b0g[1:nr] + scal$lr - dcal$lr
        E[selected,selected] <- sp$lr$cov[1:nr,1:nr]
        if (exterr){
            J[selected,ns*nr+(1:nr)] <- diag(nr)
            J[selected,(ns+1)*nr+(1:nr)] <- -diag(nr)
        }
    }
    E[ns*nr+(1:nr),ns*nr+(1:nr)] <- scal$cov
    E[(ns+1)*nr+(1:nr),(ns+1)*nr+(1:nr)] <- dcal$cov
    calibrated$cov <- J %*% E %*% t(J)
    out$calibrated <- calibrated
    class(out) <- unique(append("calibrated",class(out)))
    out
}

calibrate_geochron <- function(dat,exterr=FALSE){
    out <- dat
    type <- datatype(dat)
    if (type=="U-Pb"){
        UPb <- fractical(dat,type="U-Pb",exterr=exterr)
        PbPb <- nofractical(dat,type="U-Pb")
        out$calibrated <- mergecal(UPb,PbPb)
    } else if (type=="U-Th-Pb"){
        UPb <- fractical(dat,type="U-Pb",exterr=exterr)
        PbPb <- nofractical(dat,type="U-Th-Pb")
        ThPb <- fractical(dat,type="Th-Pb",exterr=exterr)
        out$calibrated <- mergecal(UPb,ThPb,PbPb)
    } else {
        stop("Invalid data type supplied to calibrate function.")
    }
    class(out) <- unique(append("calibrated",class(out)))
    out
}

fractical <- function(dat,type="U-Pb",exterr=FALSE){
    scal <- do.call(dat$standard$fetch,args=list(dat=dat)) # standard calibration
    dcal <- dat$calibration[[type]] # data calibration
    if (type=='U-Pb'){
        num='Pb206'
        den='U238'
    } else if (type=='Th-Pb'){
        num='Pb208'
        den='Th232'
    }
    snames <- names(dat$samples)
    ns <- length(snames)
    out <- list()
    out$num <- num
    out$den <- den
    fit <- dcal$fit
    fitcov <- matrix(c(fit$a[2]^2,fit$cov.ab,
                       fit$cov.ab,fit$b[2]^2),2,2)
    DP <- paste0(num,den)
    E <- matrix(0,ns+3,ns+3)
    E[ns+1,ns+1] <- scal$cov[DP,DP]
    E[ns+(2:3),ns+(2:3)] <- fitcov
    J <- matrix(0,ns,ns+3)
    out$lr <- rep(0,ns)
    rownames(J) <- snames
    names(out$lr) <- snames
    for (i in 1:ns){
        sp <- spot(dat,i=i)
        b0g <- sp$lr$b0g
        bXlab <- paste0('b0[',dcal$oxide,'/',den,']')
        gXlab <- paste0('g[',dcal$oxide,'/',den,']')
        bYlab <- paste0('b0[',num,'/',den,']')
        gYlab <- paste0('g[',num,'/',den,']')
        tt <- dcal$t
        XY <- rep(0,2)
        XY[1] <- b0g[bXlab] + tt*b0g[gXlab]
        XY[2] <- b0g[bYlab] + tt*b0g[gYlab]
        Eb0g <- sp$lr$cov[c(bXlab,gXlab,bYlab,gYlab),
                          c(bXlab,gXlab,bYlab,gYlab)]
        Jb0g <- matrix(0,2,4)
        Jb0g[1,1] <- 1
        Jb0g[1,2] <- tt
        Jb0g[2,3] <- 1
        Jb0g[2,4] <- tt
        E[i+(0:1),i+(0:1)] <- Jb0g %*% Eb0g %*% t(Jb0g)
        out$lr[i] <- XY[2] + scal$lr[DP] - (fit$a[1]+fit$b[1]*XY[1])
        J[i,i] <- -fit$b[1]  # dlrdX
        J[i,i+1] <- 1        # dlrdY
        if (exterr){
            J[i,ns+1] <- 1       # dlrdscalDP
            J[i,ns+2] <- -1      # dlrda
            J[i,ns+3] <- -XY[1]  # dlrdb
        }    
    }
    out$cov <- J %*% E %*% t(J)
    out
}

nofractical <- function(dat,type="U-Pb"){
    if (type=='U-Pb'){
        num=c('Pb204','Pb207')
        den=c('Pb206','Pb206')
    } else if (type=='U-Th-Pb'){
        num=c('Pb204','Pb207','Pb204')
        den=c('Pb206','Pb206','Pb208')
    }    
    snames <- names(dat$samples)
    ns <- length(snames)
    nr <- length(num)
    out <- list()
    out$num <- num
    out$den <- den
    out$lr <- rep(0,nr*ns)
    out$cov <- matrix(0,nr*ns,nr*ns)
    for (i in 1:ns){
        sp <- spot(dat,i=i)
        j <- (i-1)*nr
        lr <- paste0('b0[',num,'/',den,']')
        out$lr[j+(1:nr)] <- sp$lr$b0g[lr]
        out$cov[j+(1:nr),j+(1:nr)] <- sp$lr$cov[lr,lr]
    }
    out
}

mergecal <- function(...){
    cals <- list(...)
    num <- NULL
    den <- NULL
    for (cal in cals){
        num <- c(num,cal$num)
        den <- c(den,cal$den)
    }
    ni <- length(num)
    ns <- length(cal$lr)/length(cal$num)
    out <- list()
    out$num <- num
    out$den <- den
    out$lr <- rep(0,ni*ns)
    out$cov <- matrix(0,ni*ns,ni*ns)
    for (cal in cals){
        i <- which(num %in% cal$num)
        ii <- as.vector(sapply((0:(ns-1))*ni,'+',i))
        out$lr[ii] <- cal$lr
        out$cov[ii,ii] <- cal$cov
    }
    out
}

#' @title plot calibrated data
#' @description shows the calibrated data on a logratio plot.
#' @param x an object of class \code{calibrated}
#' @param type for U-Pb or U-Th-Pb data, either \code{'U-Pb'} or
#'     \code{'Th-Pb'}.
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function
#' @examples
#' data('Cameca')
#' dc <- drift(x=Cameca)
#' lr <- logratios(x=dc)
#' cal <- calibration(lr=lr,stand=standard(preset='Plesovice'))
#' cd <- calibrate(cal)
#' plot(cd,type='U-Pb')
#' @method plot calibrated
#' @export
plot.calibrated <- function(x,option=1,...){
    if (stable(x)) {
        out <- caldplot_stable(dat=x,...)
    } else {
        out <- caldplot_geochronology(dat=x,option=option,...)
    }
    invisible(out)
}

caldplot_stable <- function(dat,...){
    cal <- dat$calibration
    num <- dat$method$num
    den <- dat$method$den
    nn <- length(num)
    np <- nn*(nn-1)/2       # number of plot panels
    nc <- ceiling(sqrt(np)) # number of rows
    nr <- ceiling(np/nc)    # number of columns
    oldpar <- par(mfrow=c(nr,nc),mar=rep(3.5,4))
    ii <- 1
    snames <- names(dat$samples)
    for (i in 1:(nn-1)){
        for (j in (i+1):nn){
            xlab <- paste0('log[',num[i],'/',den[i],']')
            ylab <- paste0('log[',num[j],'/',den[j],']')
            B <- beta2york(lr=dat,num=num[c(i,j)],den=den[c(i,j)])
            xlim <- c(min(cal$lr[i]-3*sqrt(cal$cov[i,i]),B[,'X']-3*B[,'sX']),
                      max(cal$lr[i]+3*sqrt(cal$cov[i,i]),B[,'X']+3*B[,'sX']))
            ylim <- c(min(cal$lr[j]-3*sqrt(cal$cov[j,j]),B[,'Y']-3*B[,'sY']),
                      max(cal$lr[j]+3*sqrt(cal$cov[j,j]),B[,'Y']+3*B[,'sY']))
            plot(xlim,ylim,type='n',ann=FALSE)
            deltagrid(dat,i,j)
            IsoplotR::scatterplot(B,xlim=xlim,ylim=ylim,add=TRUE,...)
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

deltagrid <- function(dat,i,j){
    usr <- par('usr')
    xlim <- usr[1:2]
    ylim <- usr[3:4]
    xs <- dat$calibration$lr[i]
    ys <- dat$calibration$lr[j]
    di <- dat$standard$val[i]
    dj <- dat$standard$val[j]
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
    matlines(rbind(xticks,xticks),matrix(rep(ylim,nxt),ncol=nxt),lty=3,col='black')
    matlines(matrix(rep(xlim,nyt),ncol=nyt),rbind(yticks,yticks),lty=3,col='black')
    axis(side=3,at=xticks,labels=dxticks)
    mtext(expression(delta*"'"),side=3,line=2)
    axis(side=4,at=yticks,labels=dyticks)
    mtext(expression(delta*"'"),side=4,line=2)
    lines(rep(xs,2),ylim,lty=2,col='red')
    text(xs,ylim[1],labels=dat$standard$prefix,pos=4,srt=90,offset=0)
    lines(xlim,rep(ys,2),lty=2,col='red')
    text(xlim[1],ys,labels=dat$standard$prefix,pos=4,offset=0)
}

caldplot_geochronology <- function(dat,option=1,...){
    cal <- dat$calibration
    np <- length(cal)       # number of plot panels
    nc <- ceiling(sqrt(np)) # number of columns
    nr <- ceiling(np/nc)    # number of rows
    oldpar <- par(mfrow=c(nr,nc))
    ii <- 1
    for (i1 in 1:(nr-1)){
        for (i2 in (i1+1):nc){
            if (ii>np) break
            caldplot_geochronology_helper(dat,option=option,type=ii,...)
            ii <- ii + 1
        }
    }
}

caldplot_geochronology_helper <- function(dat,option=1,type=1,...){
    cal <- dat$calibration[[type]]
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
    plot(xlim,ylim,type='n',ann=FALSE)
    agegrid(dat,type,xlim,ylim)
    if (option==1){
        IsoplotR::scatterplot(yd,fit=cal$fit,xlab=xlab,ylab=ylab,
                              xlim=xlim,ylim=ylim,add=TRUE,...)
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
        IsoplotR::scatterplot(yd,fit=cal$fit,xlab=xlab,ylab=ylab,
                              xlim=xlim,ylim=ylim,add=TRUE,...)
        graphics::matlines(X,Y,lty=1,col='darkgrey')
        if (option>2){
            graphics::points(X[1,],Y[1,],pch=21,bg='black')
            graphics::points(X[nrow(X),],Y[nrow(Y),],pch=21,bg='white')
        }
    }
}

agegrid <- function(dat,type,xlim,ylim){
    cal <- dat$calibration[[type]]
    st <- do.call(dat$standard$fetchfun,args=list(dat=dat))
    lrlim <- rep(0,2)
    a <- cal$fit$a[1]
    b <- cal$fit$b[1]
    adj1 <- st$lr[type] - a - b*xlim[2]
    adj2 <- st$lr[type] - a - b*xlim[1]
    lrlim[1] <- ylim[1] + adj1
    lrlim[2] <- ylim[2] + adj2
    tlim <- IsoplotR::age(exp(lrlim),method=chronometer(dat,type))
    tticks <- pretty(tlim)
    nt <- length(tticks)
    ratio <- names(st$lr)[type]
    lrmin <- log(IsoplotR::age2ratio(tticks,ratio=ratio)[,1]) - adj2
    lrmax <- lrmin + cal$fit$b[1]*diff(xlim)
    matlines(matrix(rep(xlim,nt),nrow=2),rbind(lrmin,lrmax),lty=2,col='gray50')
}
