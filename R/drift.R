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
        ni <- length(j)
        fit <- optim(par=0,f=SS_g,method='BFGS',spot=spot,ions=ions[j])
        out['g',ions[j]] <- fit$par
        out['exp_a0',ions[j]] <- get_exp_a0(g=fit$par,spot=spot,ions=ions[j])
    }
    out
}

get_exp_a0 <- function(g,spot,ions){
    out <- rep(0,length(ions))
    names(out) <- ions
    for (ion in ions){
        p <- alphapars(spot,ion)
        if (all(p$sig>p$bkg)){
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
            out$counts <- spot$signal[,ion]*dwelltime
            out$edt <- dwelltime-deadtime*out$counts*1e-9
            out$sig <- out$counts/out$edt
            out$bkgcounts <- out$bkg*dwelltime
        } else {
            dwelltime <- spot$dwelltime[ion]
            out$counts <- spot$signal[,ion]*dwelltime
            out$edt <- dwelltime-spot$deadtime*out$counts*1e-9
            out$sig <- out$counts/out$edt
            out$bkgcounts <- out$bkg*dwelltime
        }
    }
    out
}

SS_g <- function(par,spot,ions){
    nions <- length(ions)
    g <- par
    exp_a0 <- get_exp_a0(g=g,spot=spot,ions=ions)
    tt <- hours(spot$time[,ions,drop=FALSE])
    nt <- nrow(tt)
    gt <- sweep(tt,2,g,'*')
    exp_a <- sweep(exp(gt),2,exp_a0,'*')
    bkg <- background(spot,ions)
    if (spot$m$nominalblank){
        predsig <- sweep(exp_a,2,bkg,'+')
    } else {
        predsig <- sweep(exp_a,1,bkg,'+')
    }
    D <- predsig - spot$signal[,ions]
    sum(D^2)
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
        ylab <- paste0(ion,'- b')
        graphics::plot(ap$t,sb,type='p',xlab='',ylab='',...)
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ylab,line=2)
        exp_a0 <- spot$dc['exp_a0',ion]
        g <- spot$dc['g',ion]
        predsig <- exp_a0*exp(g*ap$t)
        graphics::lines(ap$t,predsig)
    }
    graphics::par(oldpar)
}
