plot_timeresolved <- function(samp,fit=FALSE,c64=NULL){
    ions <- names(samp$dwelltime)
    np <- length(ions)      # number of plot panels
    nr <- ceiling(sqrt(np)) # number of rows
    nc <- ceiling(np/nr)    # number of columns
    graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    simplex <- c('Pb204','Pb206','Pb207','U238','UO2')
    if (fit){
        X <- samp$time[,simplex]
        Y <- predict_counts(samp,c64=c64)[,simplex]
    }
    for (ion in ions){
        graphics::plot(samp$time[,ion],samp$counts[,ion],type='p',xlab='',ylab='')
        if (fit & ion%in%simplex){
            graphics::lines(X[,ion],Y[,ion])
        }
        graphics::mtext(side=1,text='t',line=2)
        graphics::mtext(side=2,text=ion,line=2)
    }
}

predict_counts <- function(samp,c64=NULL){
    simplex <- c('Pb204','Pb206','Pb207','U238','UO2')
    Lm <- raw_count_ratios(samp=samp)
    aLm <- avg_Lm(Lm)
    ag <- get_ag(samp=samp,c64=c64)
    if (is.null(c64)){
        L4m <- log(exp(ag$x['a4'])+1)
    } else {
        d4 <- samp$dwelltime['Pb204']
        d6 <- samp$dwelltime['Pb206']
        L4m <- log(exp(ag$x['a4'])+1) + log(c64*d6/d4)
    }
    L6m <- aLm$x['L6m']
    L7m <- aLm$x['L7m']
    LUm <- aLm$x['LUm']
    den <- exp(L6m-L4m)+exp(L6m)+exp(L7m+L6m)+exp(LUm)+1
    p4 <- exp(L6m-L4m)/den
    p6 <- exp(L6m)/den
    p7 <- exp(L7m+L6m)/den
    pU <- 1/den
    pUO <- exp(LUm)/den
    p <- c(p4,p6,p7,pU,pUO)
    rs <- rowSums(samp$counts[,simplex])
    dt <- samp$dwelltime[simplex]
    counts <- matrix(rs,ncol=1) %*% matrix(p,nrow=1)
    #cps <- sweep(counts,MARGIN=2,FUN='/',dt)
    colnames(counts) <- simplex
    counts
}

calplot_raw <- function(dat,i=NULL,c64=NULL,...){
    vars <- reshuffle(dat)
    nt <- nrow(dat[[1]]$time)
    ns <- length(dat)
    if (is.null(i)) ii <- 1:ns
    else ii <- i
    graphics::matplot(vars$X[,ii],vars$Y[,ii],type='l',...)
    graphics::text(vars$X[,ii],vars$Y[,ii],labels=1:nt)
}

reshuffle <- function(samps,c64=NULL){
    out <- list(counts=list(),cps=list(),time=list(),
                sbm=list(),dwelltime=NULL)
    snames <- names(samps)
    ns <- length(snames)
    ions <- names(samps[[1]]$dwelltime)
    for (i in 1:ns){
        samp <- samps[[i]]
        for (ion in ions){
            out$counts[[ion]] <- cbind(out$counts[[ion]], samp$counts[,ion])
            out$cps[[ion]] <- cbind(out$cps[[ion]],
                                    samp$cps[,ion] + 0.5/samp$dwelltime[ion])
            out$time[[ion]] <- cbind(out$time[[ion]], samp$time[,ion])
            out$sbm[[ion]] <- cbind(out$sbm[[ion]], samp$sbm[,ion])
        }
        out$dwelltime <- rbind(out$dwelltime, samp$dwelltime)
    }
    for (ion in ions){
        colnames(out$counts[[ion]]) <- snames
        colnames(out$cps[[ion]]) <- snames
        colnames(out$time[[ion]]) <- snames
    }
    out$X <- log(out$cps[['UO2']]) - log(out$cps[['U238']])
    if (is.null(c64)){
        out$Y <- log(out$cps[['Pb206']])-log(out$cps[['U238']])
    } else {
        out$Y <- log(out$cps[['Pb206']]- out$cps[['Pb204']]*c64) -
            log(out$cps[['U238']])
    }
    colnames(out$X) <- snames
    colnames(out$Y) <- snames
    rownames(out$dwelltime) <- snames
    out
}

calplot <- function(XY,fit,alpha=0.05,sigdig=2,omit=NULL,...){
    xy <- flatXYtable(XY)
    isofit <- list(model=1)
    isofit$fact <- stats::qt(1-alpha/2,fit$df)
    isofit$n <- nrow(xy)
    isofit$mswd <- fit$mswd
    isofit$p.value <- as.numeric(1-stats::pchisq(fit$mswd/fit$df,fit$df))
    A <- fit$AB['A']
    B <- fit$AB['B']
    sA <- sqrt(fit$cov['A','A'])
    sB <- sqrt(fit$cov['B','B'])
    isofit$a <- c(A,sA,isofit$fact*sA,sqrt(isofit$mswd)*isofit$fact*sA)
    isofit$b <- c(B,sB,isofit$fact*sB,sqrt(isofit$mswd)*isofit$fact*sB)
    isofit$cov.ab <- fit$cov['A','B']
    xlab <- quote(''^238*'U'^16*'O'[2]*'/'^238*'U')
    ylab <- quote(''^206*'Pb*/'^238*'U')
    IsoplotR::scatterplot(xy,xlab=xlab,ylab=ylab,fit=isofit,omit=omit,...)
    graphics::title(IsoplotR:::isochrontitle(isofit,sigdig=sigdig),
                    xlab=xlab,ylab=ylab)
}

flatXYtable <- function(XY){
    snames <- names(XY)
    ns <- length(snames)
    xy <- matrix(0,ns,5)
    colnames(xy) <- c('X','sX','Y','sY','rXY')
    for (i in 1:ns){
        xy[i,'X'] <- XY[[i]]$x['X']
        xy[i,'Y'] <- XY[[i]]$x['Y']
        xy[i,'sX'] <- sqrt(XY[[i]]$cov['X','X'])
        xy[i,'sY'] <- sqrt(XY[[i]]$cov['Y','Y'])
        xy[i,'rXY'] <- XY[[i]]$cov['X','Y']/(xy[i,'sX']*xy[i,'sY'])
    }
    xy
}

plot_raw_ratios <- function(dat,num='Pb206',den='UO',plot=TRUE){
    snames <- names(dat)
    ns <- length(snames)
    nr <- nrow(dat[[1]]$counts)
    out <- matrix(0,nr,ns)
    colnames(out) <- snames
    for (sname in snames){
        out[,sname] <- log(dat[[sname]]$counts[,num]) -
            log(dat[[sname]]$counts[,den])
    }
    X <- matrix(1:ns,nr,ns,byrow=FALSE)
    if (plot) plot(X,out)
    else return(out)
}
plot_avg_ratios <- function(dat,num='Pb206',den='UO',radial=FALSE){
    R <- plot_raw_ratios(dat,num=num,den=den,plot=FALSE)
    avg <- colMeans(R)
    se <- apply(R,2,sd)/sqrt(nrow(R))
    XsX <- cbind(avg,se)
    if (radial) IsoplotR:::radial_helper(XsX,transformation='linear')
    else IsoplotR::weightedmean(XsX)
}

plot_raw_scatter <- function(dat,num='Pb206',den='UO'){
    snames <- names(dat)
    np <- length(snames)
    nt <- nrow(dat[[1]]$counts)
    nr <- ceiling(sqrt(np))
    nc <- ceiling(np/nr)
    dev.new(width=2*nc,height=2*nr,units = "in")
    graphics::par(mfrow=c(nr,nc),mar=c(3.5,3.5,0.5,0.5))
    for (sname in snames){
        X <- log(dat[[sname]]$counts[,num])
        Y <- log(dat[[sname]]$counts[,den])
        plot(X,Y,type='p')
        mtext(sname,line=-1)
    }
}
