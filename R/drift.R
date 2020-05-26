drift <- function(x){
    out <- x
    snames <- names(x$x)
    for (sname in snames){
        sp <- spot(dat=x,sname=sname)
        out$x[[sname]]$dc <- drift.spot(spot=sp)
    }
    class(out) <- append('drift',class(out))
    out
}

drift.spot <- function(spot){
    ions <- spot$m$ions
    nions <- length(ions)
    el <- element(ions)
    EL <- unique(el)
    nEL <- length(EL)
    out <- matrix(0,2,nions)
    colnames(out) <- ions
    rownames(out) <- c('exp_a0','g')
    for (i in 1:nEL){ # loop through the elements
        j <- which(el %in% EL[i])
        fit <- optim(par=0,fn=misfit_g,method='BFGS',spot=spot,ions=ions[j])
        out['g',ions[j]] <- fit$par
        out['exp_a0',ions[j]] <- get_exp_a0(g=fit$par,spot=spot,ions=ions[j])
    }
    out
}

misfit_g <- function(par,spot,ions){
    g <- par
    out <- 0
    exp_a0 <- get_exp_a0(g=g,spot=spot,ions=ions)
    for (ion in ions){
        p <- alphapars(spot=spot,ion=ion)
        exp_a <- exp_a0[ion]*exp(g*p$t)
        if (spot$dtype[ion]=='Fc'){
            D <-  log(exp_a) - log(p$sig - p$bkg)
            out <- out + sum(D^2)
        } else {
            obs <- p$counts
            pred <- p$bkgcounts + exp_a*p$edt
            out <- out - sum(dpois(x=obs,lambda=pred,log=TRUE))
        }
    }
    out
}

get_exp_a0 <- function(g,spot,ions){
    out <- rep(0,length(ions))
    names(out) <- ions
    for (ion in ions){
        p <- alphapars(spot,ion)
        if (spot$dtype[ion]=='Fc'){
            out[ion] <- sum(exp(g*p$t)*(p$sig-p$bkg))/sum(exp(g*p$t)^2)
        } else {
            out[ion] <- 1e-10
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
        if (spot$m$instrument=='Cameca'){
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
            out$bkgcounts <- out$bkg
        }
    }
    out
}

plot.drift <- function(x,sname,i=1,...){
    spot <- spot(x,sname=sname,i=i)
    ions <- spot$m$ions
    np <- length(ions) # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    oldpar <- graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (ion in ions){
        ap <- alphapars(spot,ion)
        sb <- ap$sig - ap$bkg
        tt <- seconds(ap$t)
        ylab <- paste0(ion,'- b')
        graphics::plot(tt,sb,type='p',xlab='',ylab='',...)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
        exp_a0 <- spot$dc['exp_a0',ion]
        g <- spot$dc['g',ion]
        predsig <- exp_a0*exp(g*ap$t)
        graphics::lines(tt,predsig)
    }
    graphics::par(oldpar)
}
