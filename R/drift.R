#' @title drift correction
#' @description fits a generalised linear model to time resolved SIMS data
#' @param x an object of class \code{simplex}
#' @return an object of class \code{drift}
#' @examples
#' data('SHRIMP',package='simplex')
#' dc <- drift(x=SHRIMP)
#' plot(dc,i=1)
#' @export
drift <- function(x){
    out <- x
    snames <- names(x$samples)
    for (sname in snames){
        sp <- spot(dat=x,sname=sname)
        out$samples[[sname]]$dc <- drift.spot(spot=sp)
    }
    class(out) <- append('drift',class(out))
    out
}

drift.spot <- function(spot){
    ions <- spot$method$ions
    nions <- length(ions)
    el <- element(ions)
    EL <- unique(el)
    nEL <- length(EL)
    out <- matrix(0,2,nions)
    colnames(out) <- ions
    rownames(out) <- c('a0','g')
    for (i in 1:nEL){ # loop through the elements
        j <- which(el %in% EL[i])
        fit <- stats::optim(par=0,fn=misfit_g,method='L-BFGS-B',
                            lower=-5,upper=5,spot=spot,ions=ions[j])
        out['g',ions[j]] <- fit$par
        out['a0',ions[j]] <- get_a0(g=fit$par,spot=spot,ions=ions[j])
    }
    out
}

misfit_g <- function(par,spot,ions){
    g <- par
    out <- 0
    a0 <- get_a0(g=g,spot=spot,ions=ions)
    for (ion in ions){
        p <- alphapars(spot=spot,ion=ion)
        a <- exp(a0[ion]+g*p$t)
        if (spot$dtype[ion]=='Fc'){
            D <-  a - log(p$sig - p$bkg)
            out <- out + sum(D^2)
        } else {
            obs <- p$counts
            pred <- a*p$edt # + p$bkgcounts # approximate
            out <- out - sum(stats::dpois(x=obs,lambda=pred,log=TRUE))
        }
    }
    out
}

get_a0 <- function(g,spot,ions){
    sem_misfit <- function(par,g,ap){
        expected <- ap$bkgcounts + exp(par+g*ap$t)*ap$edt
        LL <- dpois(x=ap$counts,lambda=expected,log=TRUE)
        -sum(LL)
    }
    out <- rep(0,length(ions))
    names(out) <- ions
    for (ion in ions){
        ap <- alphapars(spot,ion)
        if (spot$dtype[ion]=='Fc'){
            out[ion] <- log(sum(exp(g*ap$t)*(ap$sig-ap$bkg))) -
                log(sum(exp(g*ap$t)^2))
        } else {
            init <- log(max(0.5,sum(ap$counts))) -
                log(sum(exp(g*ap$t)*ap$edt))
            out[ion] <- optimise(sem_misfit,interval=init+c(-5,2),
                                 g=g,ap=ap)$minimum
        }
    }
    out
}

alphapars <- function(spot,ion){
    out <- list()
    out$bkg <- background(spot,ion)
    out$t <- hours(spot$time[,ion])
    if (faraday(spot,ion)){
        out$sig <- spot$signal[,ion]
    } else {
        if (spot$method$instrument=='Cameca'){
            detector <- spot$detector[ion]
            deadtime <- spot$deadtime[detector]
            dwelltime <- spot$dwelltime[ion]
            counts <- spot$signal[,ion]*dwelltime
            out$edt <- dwelltime-deadtime*counts*1e-9
            out$sig <- counts/out$edt
            out$counts <- round(counts)
            out$bkgcounts <- round(out$bkg*dwelltime)
        } else {
            dwelltime <- spot$dwelltime[ion]
            out$counts <- spot$signal[,ion]
            out$edt <- dwelltime-spot$deadtime*out$counts*1e-9
            out$sig <- out$counts/out$edt
            out$bkgcounts <- spot$signal[,'bkg']
        }
    }
    out
}

#' @title plot drift corrected data
#' @description shows the time resolved mass spectrometer signals
#'     fitted by a generalised linear model
#' @param x an object of class \code{drift}
#' @param sname the name of the sample to be shown
#' @param i the number of the sample to be shown
#' @param ... optional arguments to be passed on to the generic
#'     \code{plot} function.
#' @examples
#' data('Cameca',package='simplex')
#' dc <- drift(x=Cameca)
#' plot(dc,i=1)
#' @method plot drift
#' @export
plot.drift <- function(x,sname=NULL,i=1,...){
    spot <- spot(x,sname=sname,i=i)
    ions <- spot$method$ions
    np <- length(ions) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in ions){
        ap <- alphapars(spot,ion)
        sb <- ap$sig - ap$bkg
        tt <- seconds(ap$t)
        ylab <- paste0(ion,'- b')
        a0 <- spot$dc['a0',ion]
        g <- spot$dc['g',ion]
        predsig <- exp(a0+g*ap$t)
        graphics::plot(rep(tt,2),c(sb,predsig),
                       type='n',xlab='',ylab='',...)
        graphics::points(tt,sb,type='p')
        graphics::lines(tt,predsig)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
    }
    graphics::par(oldpar)
}
