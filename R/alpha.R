alpha <- function(spot,ions=spot$ions,plot=FALSE){
    init_a0 <- function(spot,ions){
        sb1 <- subtract_blank(spot=spot,ions=ions)[1,]
        sb1[sb1<0] <- -sb1[sb1<0]
        log(sb1)
    }
    nions <- length(ions)
    el <- element(ions)
    EL <- unique(el)
    nEL <- length(EL)
    ia0 <- init_a0(spot=spot,ions=ions)
    a0 <- rep(0,nions)
    names(a0) <- ions
    g <- rep(0,nEL)
    names(g) <- EL
    for (i in 1:nEL){ # loop through elements
        j <- which(el %in% EL[i])
        gfit <- optim(par=0,f=LL_g,method='L-BFGS-B',lower=-1,upper=1,
                      spot=spot,ions=ions[j],ia0=ia0[ions[j]])
        afit <- optim(par=ia0,fn=LL_a0,method='BFGS',gr=NULL,
                      spot=spot,ions=ions[j],g=gfit$par)
        g[EL[i]] <- gfit$par
        a0[ions[j]] <- afit$par
    }
    out <- list(a0=a0,g=g)
    if (plot){
        plot_alpha(spot=spot,a0g=out)
    }
    out
}

LL_g <- function(par,spot,ions=spot$ions,ia0){
    fit <- optim(par=ia0,fn=LL_a0,gr=NULL,spot=spot,ions=ions,g=par)
    print(fit$value)
#    plot_alpha(spot=spot,ions=ions,a0g=list(a0=fit$par,g=par))
    fit$value
}

LL_a0 <- function(par,spot,ions=spot$ions,g){
    LL_a0g(c(par,g),spot=spot,ions=ions)
}

LL_a0g <- function(a0g,spot,ions=spot$ions){
    nions <- length(ions)
    a0 <- a0g[1:nions]
    g <- a0g[nions+1]
    tt <- hours(spot$time[,ions,drop=FALSE])
    nt <- nrow(tt)
    gt <- sweep(tt,2,g,'*')
    a <- sweep(gt,2,a0,'+')
    detector <- spot$detector[ions]
    if (spot$nominalblank){
        bkg <- spot$background[detector]
        predsig <- sweep(exp(a),2,bkg,'+')
        E <- diag(nions)
        D <- predsig - spot$signal[,ions]
#        ij <- expand.grid(1:nions,1:nions)
#        for (r in 1:(nions^2)){
#            i <- ij[r,1]
#            j <- ij[r,2]
#            E[i,j] <- sum(D[,i]*D[,j])/(nt-1)
#        }
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
