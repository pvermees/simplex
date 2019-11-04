logratios <- function(dat,c64=NULL){
    snames <- names(dat)
    ns <- length(snames)
    out <- list()
    for (sname in snames){
        print(sname)
        samp <- dat[[sname]]
        Lm <- raw_count_ratios(samp=samp)
        aLm <- avg_Lm(Lm)
        ag <- get_ag(samp=samp,c64=c64)
        out[[sname]] <- bias_correction(samp=samp,aLm=aLm,ag=ag)
    }
    out
}

raw_count_ratios <- function(samp){
    counts <- samp$counts[,c('206Pb','207Pb','238U','238U16O2')]
    lc <- log(counts)
    nt <- nrow(lc)
    n6 <- counts[,'206Pb']
    n7 <- counts[,'207Pb']
    nU <- counts[,'238U']
    nUO <- counts[,'238U16O2']
    L6m <- lc[,'206Pb'] - lc[,'238U']
    L7m <- lc[,'207Pb'] - lc[,'206Pb']
    LUm <- lc[,'238U16O2'] - lc[,'238U']
    rs <- rowSums(counts)
    D <- exp(L6m)+exp(L7m+L6m)+exp(LUm)+1
    dDdL6m <- exp(L6m)+exp(L7m+L6m)
    dDdL7m <- exp(L7m+L6m)
    dDdLUm <- exp(LUm)
    d2DdL6m2 <- exp(L6m)+exp(L7m+L6m)
    d2DdL7m2 <- exp(L7m+L6m)
    d2DdLUm2 <- exp(LUm)
    d2DdL6mdL7m <- exp(L7m+L6m)
    d2DdL7mdL6m <- exp(L7m+L6m)
    d2DdL6mdLUm <- 0
    d2DdL7mdLUm <- 0
    d2DdLUmdL6m <- 0
    d2DdLUmdL7m <- 0
    #LL <- n6*L6m + n7*(L6m+L7m) + nUO*LUm - rs*log(D)
    d2LLdL6m2 <- (rs/D^2)*(dDdL6m)^2 - (rs/D)*d2DdL6m2
    d2LLdL6mdL7m <- (rs/D^2)*dDdL6m*dDdL7m - (rs/D)*d2DdL6mdL7m
    d2LLdL6mdLUm <- (rs/D^2)*dDdL6m*dDdLUm - (rs/D)*d2DdL6mdLUm
    d2LLdL7m2 <- (rs/D^2)*(dDdL7m)^2 - (rs/D)*d2DdL7m2
    d2LLdL7mdL6m <- (rs/D^2)*dDdL7m*dDdL6m - (rs/D)*d2DdL6mdL7m
    d2LLdL7mdLUm <- (rs/D^2)*dDdL7m*dDdLUm - (rs/D)*d2DdL7mdLUm
    d2LLdLUm2 <- (rs/D^2)*(dDdLUm)^2 - (rs/D)*d2DdLUm2
    d2LLdLUmdL6m <- (rs/D^2)*dDdLUm*dDdL6m - (rs/D)*d2DdLUmdL6m
    d2LLdLUmdL7m <- (rs/D^2)*dDdLUm*dDdL7m - (rs/D)*d2DdLUmdL7m
    out <- list()
    labels <- c('L6m','L7m','LUm')
    for (i in 1:nt){
        out[[i]] <- list()
        out[[i]]$x <- c(L6m[i],L7m[i],LUm[i])
        H <- matrix(0,3,3)
        H[1,1] <- d2LLdL6m2[i]
        H[1,2] <- d2LLdL6mdL7m[i]
        H[1,3] <- d2LLdL6mdLUm[i]
        H[2,1] <- d2LLdL7mdL6m[i]
        H[2,2] <- d2LLdL7m2[i]
        H[2,3] <- d2LLdL7mdLUm[i]
        H[3,1] <- d2LLdLUmdL6m[i]
        H[3,2] <- d2LLdLUmdL7m[i]
        H[3,3] <- d2LLdLUm2[i]
        out[[i]]$H <- H
        names(out[[i]]$x) <- labels
        colnames(out[[i]]$H) <- labels
        rownames(out[[i]]$H) <- labels
    }
    out
}

avg_Lm <- function(Lm){
    np <- length(Lm[[1]]$x)
    init <- rep(0,np)
    nt <- length(Lm)
    for (i in 1:nt){ init <- init + Lm[[i]]$x/nt }
    lower <- init - rep(1,np)
    upper <- init + rep(1,np)
    fit <- stats::optim(init,misfit_avg_Lm,Lm=Lm,
                        lower=lower,upper=upper,
                        method='L-BFGS-B',hessian=TRUE)
    out <- list()
    out$x <- fit$par
    out$cov <- solve(fit$hessian)
    labels <- c('L6m','L7m','LUm')
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}
misfit_avg_Lm <- function(aLm,Lm){
    nt <- length(Lm)
    out <- 0
    for (i in 1:nt){
        D <- Lm[[i]]$x - aLm
        out <- out - D %*% Lm[[i]]$H %*% D
    }
    out/2
}

get_ag <- function(samp,c64=NULL){
    init <- init_ag(samp,c64=c64)
    np <- length(init) # number of parameters
    lower <- init - rep(1,np)
    upper <- init + rep(1,np)
    fit <- stats::optim(init,misfit_ag,method='L-BFGS-B',
                        lower=lower,upper=upper,samp=samp,
                        c64=c64,control=list(fnscale=-1),hessian=TRUE)
    out <- list()
    out$x <- fit$par
    out$cov <- solve(-fit$hessian)
    out    
}
init_ag <- function(samp,c64=NULL){
    t4 <- hours(samp$time[,'204Pb'])
    t6 <- hours(samp$time[,'206Pb'])
    c4 <- samp$counts[,'204Pb']
    c6 <- samp$counts[,'206Pb']
    fit6 <- stats::glm(c6 ~ t6, family=stats::poisson(link='log'))
    b6 <- fit6$coef[1]
    g <- fit6$coef[2]
    b4 <- stats::glm(c4 ~ 1 + offset(g*t4),
                     family=stats::poisson(link='log'))$coef
    if (is.null(c64)){
        a4 <- log(exp(b6-b4)-1)
    } else {
        d4 <- samp$dwelltime['204Pb']
        d6 <- samp$dwelltime['206Pb']
        a4 <- log(exp(b6-b4)*(d4/d6)/c64 - 1)
    }
    out <- c(a4,g)
    names(out) <- c('a4','g')
    out
}
misfit_ag <- function(ag,samp,c64=NULL){
    simplex <- c('204Pb','206Pb','207Pb','238U','238U16O2')
    b <- ag2b(ag,samp=samp,c64=c64)
    p <- b2p(b,ag['g'],tt=hours(samp$time[,simplex]))
    counts <- samp$counts[,simplex]
    stats::dmultinom(counts,prob=p,log=TRUE)
}
ag2b <- function(ag,samp,c64=NULL){
    simplex <- c('204Pb','206Pb','207Pb','238U','238U16O2')
    tt <- hours(samp$time[,simplex])
    cc <- samp$counts[,simplex]
    dt <- samp$dwelltime[simplex]
    fn6 <- cc[,'206Pb'] ~ 1 + offset(ag['g']*tt[,'206Pb'])
    fn7 <- cc[,'207Pb'] ~ 1 + offset(ag['g']*tt[,'207Pb'])
    fnU <- cc[,'238U'] ~ 1 + offset(ag['g']*tt[,'238U'])
    fnUO <- cc[,'238U16O2'] ~ 1 + offset(ag['g']*tt[,'238U16O2'])
    b6 <- stats::glm(fn6,family=stats::poisson(link='log'))$coef
    b7 <- stats::glm(fn7,family=stats::poisson(link='log'))$coef
    bU <- stats::glm(fnU,family=stats::poisson(link='log'))$coef
    bUO <- stats::glm(fnUO,family=stats::poisson(link='log'))$coef
    if (is.null(c64)){
        b4 <- b6 - log(exp(ag['a4'])+1)
    } else {
        b4 <- b6 - log(c64*dt['206Pb']/dt['204Pb']) - log(exp(ag['a4'])+1)
    }
    out <- c(b4,b6,b7,bU,bUO)
    names(out) <- simplex
    out
}
b2p <- function(b,ag,tt){
    simplex <- c('204Pb','206Pb','207Pb','238U','238U16O2')
    nr <- nrow(tt)
    bm <- matrix(rep(b,nr),nrow=nr,byrow=TRUE)
    colnames(bm) <- simplex
    ebgt <- exp(bm + ag['g']*tt)
    out <- ebgt/sum(ebgt)
    colnames(out) <- simplex
    out
}

bias_correction <- function(samp,aLm,ag,c64=NULL){
    simplex <- c('204Pb','206Pb','207Pb','238U','238U16O2')
    tt <- hours(samp$time[,simplex])
    dt6 <- mean(tt[,'206Pb']-tt[,'238U'])
    dt7 <- mean(tt[,'207Pb']-tt[,'206Pb'])
    dtU <- mean(tt[,'238U16O2']-tt[,'238U'])
    L6c <- aLm$x['L6m'] - ag$x['g']*dt6
    L7c <- aLm$x['L7m'] - ag$x['g']*dt7
    LUc <- aLm$x['LUm'] - ag$x['g']*dtU
    if (is.null(c64)){
        L4c <- log(exp(ag$x['a4'])+1)
    } else {
        d4 <- samp$dwelltime['204Pb']
        d6 <- samp$dwelltime['206Pb']
        L4c <- log(exp(ag$x['a4'])+1) + log(c64*d6/d4)
    }
    out <- list()
    out$x <- c(L4c,L6c,L7c,LUc)
    labels <- c('L4c','L6c','L7c','LUc')
    names(out$x) <- labels
    # joint covariance matrix of aLm and ag:
    E <- matrix(0,5,5)
    E[1:3,1:3] <- aLm$cov
    E[4:5,4:5] <- ag$cov
    J <- matrix(0,4,5)
    rownames(J) <- labels
    colnames(J) <- c('L6m','L7m','LUm','a4','g')
    J['L6c','L6m'] <- 1
    J['L7c','L7m'] <- 1
    J['LUc','LUm'] <- 1
    J['L4c','a4'] <- exp(ag$x['a4'])/(exp(ag$x['a4'])+1) # dL4c/da4
    J['L6c','g'] <- -dt6
    J['L7c','g'] <- -dt7
    J['LUc','g'] <- -dtU
    out$cov <- J %*% E %*% t(J)
    out
}

lr2XY <- function(lr,dat,c64=NULL){
    snames <- names(dat)
    ns <- length(snames)
    out <- list()
    labels <- c('X','Y','L7a','L4a')
    calions <- c('L4c','L6c','L7c','LUc')
    J <- matrix(0,4,4)
    rownames(J) <- labels
    colnames(J) <- calions
    J['X','LUc'] <- 1
    J['Y','L6c'] <- 1
    for (sname in snames){
        d4 <- dat[[sname]]$dwelltime['204Pb']
        d6 <- dat[[sname]]$dwelltime['206Pb']
        d7 <- dat[[sname]]$dwelltime['207Pb']
        dU <- dat[[sname]]$dwelltime['238U']
        dUO <- dat[[sname]]$dwelltime['238U16O2']
        L4c <- lr[[sname]]$x['L4c']
        L6c <- lr[[sname]]$x['L6c']
        L7c <- lr[[sname]]$x['L7c']
        LUc <- lr[[sname]]$x['LUc']
        X <- LUc + log(dU/dUO)
        if (is.null(c64)){
            Y <- L6c + log(dU/d6)
        } else {
            Y <- log(dU) + log(dU/d6) + log(1-(c64*d6)/(d4*exp(L4c)))
            J['Y','L4c'] <- c64*d6/(c64*d6 + exp(L4c)*d4)
        }
        L7a <- L7c + log(d6/d7)
        L4a <- L4c + log(d6/d4)
        J['L7a','L7c'] <- log(d6/d7)
        J['L4a','L4c'] <- log(d6/d4)
        out[[sname]] <- list()
        out[[sname]]$x <- c(X,Y,L7a,L4a)
        names(out[[sname]]$x) <- labels
        out[[sname]]$cov <- J %*% lr[[sname]]$cov[calions,calions] %*% t(J)
        out[[sname]]$omega <- solve(out[[sname]]$cov)
    }
    out
}

# build a calibration curve
calibration <- function(lr,dat,plot=TRUE,disp=TRUE,omit=NULL,...){
    XY <- lr2XY(lr=lr,dat=dat)
    fit <- yorkfit(XY,omit=omit)
    out <- list()
    out$AB <- fit$par
    out$df <- length(dat)-2
    out$mswd <- 2*fit$value/out$df
    if (disp) d <- out$mswd
    else d <- 1
    out$cov <- d*solve(fit$hessian)
    labels <- c('A','B')
    names(out$AB) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    if (plot) calplot(XY,fit=out,omit=omit)
    out
}

yorkfit <- function(XY,omit=NULL){
    if (!is.null(omit)){
        keep <- !((1:length(XY))%in%omit)
        XY <- subset(XY,subset=keep)
    }
    snames <- names(XY)
    ns <- length(snames)
    X <- rep(0,ns)
    Y <- rep(0,ns)
    for (i in 1:ns){
        X[i] <- XY[[i]]$x['X']
        Y[i] <- XY[[i]]$x['Y']
    }
    init <- stats::lm(Y ~ X)$coef
    stats::optim(init,misfit_york,method='BFGS',XY=XY,hessian=TRUE)
}

misfit_york <- function(AB,XY=XY){
    A <- AB[1]
    B <- AB[2]
    SS <- 0
    snames <- names(XY)
    D <- matrix(0,1,2)
    for (sname in snames){
        X <- XY[[sname]]$x['X']
        Y <- XY[[sname]]$x['Y']
        O <- XY[[sname]]$omega[c('X','Y'),c('X','Y')]
        CC <- Y - A - B*X
        C1 <- O[1,1] + O[1,2]*B + O[2,1]*B + O[2,2]*B^2
        C2 <- 2*(O[1,2] + O[2,2]*B)*CC
        C3 <- O[2,2]*CC^2
        K <- -C2/(2*C1)
        SS <- SS + C1*K^2 + C2*K + C3
    }
    SS/2
}

# apply a calibration to a unknown samples
calibrate <- function(lr,fit,dat,PbUstand=NULL,tst=c(337.13,0.18)){
    XY <- lr2XY(lr=lr,dat=dat)
    ns <- length(XY)
    if (is.null(PbUstand))
        Pb6U8stand <- IsoplotR:::age_to_Pb206U238_ratio(tt=tst[1],st=tst[2])
    Ytstand <- log(Pb6U8stand[1])
    sYtstand <- Pb6U8stand[2]/Pb6U8stand[1]
    La <- rep(0,3*ns)
    E <- matrix(0,4*ns+3,4*ns+3)
    E[4*ns+(1:2),4*ns+(1:2)] <- fit$cov
    E[4*ns+3,4*ns+3] <- sYtstand^2
    J <- matrix(0,nrow=3*ns,ncol=4*ns+3)
    for (i in 1:ns){
        i1 <- i      # X
        i2 <- i+ns   # Y
        i3 <- i+2*ns # L7a
        i4 <- i+3*ns # L4a
        o1 <- i      # out: L6a
        o2 <- i+ns   # out: L7a
        o3 <- i+2*ns # out: L4a
        xy <- XY[[i]]
        La[o1] <- xy$x['Y'] + Ytstand - fit$AB['A'] - fit$AB['B']*xy$x['X']
        La[o2] <- xy$x['L7a']
        La[o3] <- xy$x['L4a']
        J[o1,i1] <- -fit$AB['B']   # dL6dX
        J[o1,i2] <- 1              # dL6dY
        J[o1,4*ns+1] <- -1         # dL6dA
        J[o1,4*ns+2] <- -xy$x['X'] # dL6dB
        J[o1,4*ns+3] <- 1          # dL6dYtstand
        J[o2,i3] <- 1              # dL7dL7
        J[o3,i4] <- 1              # dL4dL4
        E[c(i1,i2,i3,i4),c(i1,i2,i3,i4)] <- xy$cov
    }
    out <- list()
    out$snames <- names(lr)
    out$labels <- c('L6a','L7a','L4a')
    out$x <- La
    out$cov <- J %*% E %*% t(J)
    out
}

logratios2ratios <- function(alr){
    ns <- length(alr$snames)
    ni <- length(alr$labels)
    out <- alr
    out$labels <- c('U8Pb6','Pb76','Pb46')
    J <- matrix(0,3*ns,3*ns)
    for (i in 1:ns){
        i1 <- i      # L6a
        i2 <- i+ns   # L7a
        i3 <- i+2*ns # L4a
        out$x[i1] <- exp(-alr$x[i1])
        out$x[i2] <- exp(alr$x[i2])
        out$x[i3] <- exp(-alr$x[i3])
        J[i1,i1] <- -exp(-alr$x[i1])
        J[i2,i2] <- exp(alr$x[i2])
        J[i3,i3] <- -exp(-alr$x[i3])
    }
    out$cov <- J %*% alr$cov %*% t(J)
    out
}

hours <- function(tt){
    tt/3600
}
