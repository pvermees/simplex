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
    rownames(out) <- c('exp_a0','g')
    for (i in 1:nEL){ # loop through the elements
        j <- which(el %in% EL[i])
        ni <- length(j)
        fit <- optim(par=0,f=SS_g,method='BFGS',spot=x,ions=ions[j])
        out['g',ions[j]] <- fit$par
        out['exp_a0',ions[j]] <- get_exp_a0(g=fit$par,spot=x,ions=ions[j])
    }
    if (plot){
        plot_alpha(spot=x,ea0g=out,ions=ions,...)
    }
    out
}

get_exp_a0 <- function(g,spot,ions){
    out <- rep(0,length(ions))
    names(out) <- ions
    for (ion in ions){
        p <- alphapars(spot,ion)
        out[ion] <- sum(exp(g*p$t)*(p$sig-p$bkg))/sum(exp(g*p$t)^2)
    }
    out
}

alphapars <- function(spot,ion){
    out <- list()
    out$sig <- spot$signal[,ion]
    out$bkg <- background(spot,ion)
    out$t <- hours(spot$time[,ion])
    out
}

SS_g <- function(par,spot,ions=spot$ions){
    nions <- length(ions)
    g <- par
    exp_a0 <- get_exp_a0(g=g,spot=spot,ions=ions)
    tt <- hours(spot$time[,ions,drop=FALSE])
    nt <- nrow(tt)
    gt <- sweep(tt,2,g,'*')
    exp_a <- sweep(exp(gt),2,exp_a0,'*')
    bkg <- background(spot,ions)
    if (spot$nominalblank){
        predsig <- sweep(exp_a,2,bkg,'+')
    } else {
        predsig <- sweep(exp_a,1,bkg,'+')
    }
    D <- predsig - spot$signal[,ions]
    sum(D^2)
}

plot_alpha <- function(spot,ions=spot$ions,ea0g,...){
    np <- length(spot$ions) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in spot$ions){
        ap <- alphapars(spot,ion)
        sb <- ap$sig - ap$bkg
        ylab <- paste0(ion,'- b')
        graphics::plot(ap$t,sb,type='p',xlab='',ylab='',...)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
        if (ion %in% ions){
            exp_a0 <- ea0g['exp_a0',ion]
            g <- ea0g['g',ion]
            predsig <- exp_a0*exp(g*ap$t)
            graphics::lines(ap$t,predsig)
        }
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
