alpha <- function(spot,ions=spot$ions,plot=FALSE){
    init_a0 <- function(spot,ions){
        sb <- subtract_blank(spot=spot,ions=ions)
        cm <- colMeans(sb)
        cm[cm<0] <- -cm[cm<0]
        log(cm)
    }
    nions <- length(ions)
    el <- element(ions)
    EL <- unique(el)
    nEL <- length(EL)
    a0 <- init_a0(spot,ions)
    g <- rep(0,nEL)
    a0g <- c(a0,g)
    gi <- rep(0,nions)
    for (i in 1:nEL){
        j <- which(el %in% EL[i])
        gi[j] <- i
    }
    fit <- optim(par=a0g,fn=LL_alpha,spot=spot,ions=ions,gi=gi)
    a0 <- fit$par[1:nions]
    names(a0) <- ions
    g <- fit$par[(nions+1):(nions+nEL)]
    names(g) <- EL
    a0g <- list(a0=a0,g=g)
    if (plot){
        plot_alpha(spot=spot,a0g=a0g)
    }
    fit$par
}

LL_alpha <- function(par,spot,ions=spot$ions,gi){
    nions <- length(ions)
    a0 <- par[1:nions]
    np <- length(par)
    g <- par[(nions+1):np][gi]
    tt <- hours(spot$time[,ions])
    nt <- nrow(tt)
    gt <- sweep(tt,2,g,'*')
    a <- sweep(gt,2,a0,'+')
    detector <- spot$detector[ions]
    if (spot$nominalblank){
        bkg <- spot$background[detector]
        predsig <- sweep(exp(a),2,bkg,'+')
        ij <- expand.grid(1:nions,1:nions)
        E <- matrix(0,nions,nions)
        D <- predsig - spot$signal[,ions]
        for (r in 1:(nions^2)){
            i <- ij[r,1]
            j <- ij[r,2]
            E[i,j] <- sum(D[,i]*D[,j])/(nt-1)
        }
        SS <- sum(D %*% solve(E) %*% t(D))
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
#    if (has_outliers(spot)){
#        nt <- nrow(spot$time)
#        pch <- rep(4,nt)
#        pch[spot$inliers] <- 1
#    } else {
        pch <- 1
#    }
    for (ion in spot$ions){
        tt <- hours(spot$time[,ion])
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
        graphics::plot(tt,sb,type='p',xlab='',ylab='',pch=pch,...)
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
