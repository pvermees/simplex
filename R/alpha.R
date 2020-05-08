#' @rdname alpha
#' @export
alpha <- function(x,...){ UseMethod("alpha",x) }
#' @rdname alpha
#' @export
alpha.default <- function(x,...){ stop('No default method.') }
#' @rdname alpha
#' @export
alpha.simplex <- function(x,ions,...){
    snames <- names(x)
    out <- list()
    for (sname in snames){
        out[[sname]] <- alpha(x=x[[sname]],ions=ions,...)
    }
    out
}
#' @rdname alpha
#' @export
alpha.standards <- function(x,ions,...){
    alpha(x=x$x,ions=ions,...)
}
#' @rdname alpha
#' @export
alpha.spot <- function(x,ions=x$ions,plot=FALSE,...){
    nions <- length(ions)
    el <- element(ions)
    EL <- unique(el)
    nEL <- length(EL)
    out <- matrix(0,2,nions)
    colnames(out) <- ions
    rownames(out) <- c('a0','g')
    for (i in 1:nEL){ # loop through the elements
        j <- which(el %in% EL[i])
        ni <- length(j)
        init <- rep(0,ni+1)
        fit <- optim(par=init,f=SS_a0g,method='BFGS',spot=x,ions=ions[j])
        out['a0',ions[j]] <- fit$par[1:ni]
        out['g',ions[j]] <- fit$par[ni+1]
    }
    if (plot){
        plot_alpha(spot=x,a0g=out,ions=ions,...)
    }
    out
}

SS_a0g <- function(a0g,spot,ions=spot$ions){
    nions <- length(ions)
    a0 <- a0g[1:nions]
    g <- a0g[nions+1]
    tt <- hours(spot$time[,ions,drop=FALSE])
    nt <- nrow(tt)
    gt <- sweep(tt,2,g,'*')
    a <- sweep(gt,2,a0,'+')
    bkg <- background(spot,ions)
    if (spot$nominalblank){
        predsig <- sweep(exp(a),2,bkg,'+')
    } else {
        predsig <- sweep(exp(a),1,bkg,'+')
    }
    D <- predsig - spot$signal[,ions]
    sum(D^2)
}

plot_alpha <- function(spot,ions=spot$ions,a0g,...){
    np <- length(spot$ions) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in spot$ions){
        tt <- hours(spot$time[,ion])
        if (ion %in% ions){
            a0 <- a0g['a0',ion]
            g <- a0g['g',ion]
            predsig <- exp(a0 + g*tt)
            bkg <- background(spot,ions)
            sweep(spot$signal[,ions,drop=FALSE],2,bkg,'-')
            ylab <- paste0('signal - blank (',ion,')')
        } else {
            sb <- spot$signal[,ion]
            ylab <- paste0('signal (',ion,')')
        }
        graphics::plot(tt,sb,type='p',xlab='',ylab='',...)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
        if (ion %in% ions) graphics::lines(tt,predsig)
    }
    graphics::par(oldpar)
}

background <- function(spot,ions){
    detector <- spot$detector[ions]
    if (spot$nominalblank){
        out <- spot$background[detector]
    } else {
        out <- spot$signal[,'bkg']
    }
    out
}
