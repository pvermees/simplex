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
        out[[sname]] <- alpha_spot(spot=x[[sname]],ions=ions,...)
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
alpha_spot <- function(spot,ions=spot$ions,plot=FALSE,...){
    nions <- length(ions)
    el <- element(ions)
    EL <- unique(el)
    nEL <- length(EL)
    a0 <- rep(0,nions)
    names(a0) <- ions
    g <- rep(0,nEL)
    names(g) <- EL
    for (i in 1:nEL){ # loop through the elements
        j <- which(el %in% EL[i])
        fit <- optim(par=c(a0[ions[j]],g[i]),f=LL_a0g,
                     method='BFGS',spot=spot,ions=ions[j])
        a0[ions[j]] <- fit$par[1:length(j)]
        g[EL[i]] <- fit$par[length(j)+1]
    }
    out <- list(a0=a0,g=g)
    if (plot){
        plot_alpha(spot=spot,a0g=out,...)
    }
    out
}

LL_a0g <- function(a0g,spot,ions=spot$ions){
    nions <- length(ions)
    a0 <- a0g[1:nions]
    g <- a0g[nions+1]
    tt <- days(spot$time[,ions,drop=FALSE])
    nt <- nrow(tt)
    gt <- sweep(tt,2,g,'*')
    a <- sweep(gt,2,a0,'+')
    detector <- spot$detector[ions]
    if (spot$nominalblank){
        bkg <- spot$background[detector]
        predsig <- sweep(exp(a),2,bkg,'+')
        D <- predsig - spot$signal[,ions]
        SS <- sum(D^2)
    } else {
        stop('Not implemented yet.')
    }
    sum(SS)
}

plot_alpha <- function(spot,ions=spot$ions,a0g,...){
    np <- length(spot$ions) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in spot$ions){
        tt <- days(spot$time[,ion])
        if (ion %in% ions){
            a0 <- a0g$a0[ion]
            el <- element(ion)
            g <- a0g$g[el]
            predsig <- exp(a0 + g*tt)
            sb <- subtract_blank(spot=spot,ions=ion)
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

subtract_blank <- function(spot,ions){
    detectors <- spot$detector[ions]
    bkg <- spot$background[detectors]
    sweep(spot$signal[,ions,drop=FALSE],2,bkg,'-')
}
