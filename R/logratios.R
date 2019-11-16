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
    counts <- samp$counts[,c('Pb206','Pb207','U238','UO2')]
    lc <- log(counts)
    nt <- nrow(lc)
    n6 <- counts[,'Pb206']
    n7 <- counts[,'Pb207']
    nU <- counts[,'U238']
    nUO <- counts[,'UO2']
    L6m <- lc[,'Pb206'] - lc[,'U238']
    L7m <- lc[,'Pb207'] - lc[,'Pb206']
    LUm <- lc[,'UO2'] - lc[,'U238']
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
    t4 <- hours(samp$time[,'Pb204'])
    t6 <- hours(samp$time[,'Pb206'])
    c4 <- samp$counts[,'Pb204']
    c6 <- samp$counts[,'Pb206']
    fit6 <- stats::glm(c6 ~ t6, family=stats::poisson(link='log'))
    b6 <- fit6$coef[1]
    g <- fit6$coef[2]
    b4 <- stats::glm(c4 ~ 1 + offset(g*t4),
                     family=stats::poisson(link='log'))$coef
    if (is.null(c64)){
        a4 <- log(exp(b6-b4)-1)
    } else {
        d4 <- samp$dwelltime['Pb204']
        d6 <- samp$dwelltime['Pb206']
        a4 <- log(exp(b6-b4)*(d4/d6)/c64 - 1)
    }
    out <- c(a4,g)
    names(out) <- c('a4','g')
    out
}
misfit_ag <- function(ag,samp,c64=NULL){
    simplex <- c('Pb204','Pb206','Pb207','U238','UO2')
    b <- ag2b(ag,samp=samp,c64=c64)
    p <- b2p(b,ag['g'],tt=hours(samp$time[,simplex]))
    counts <- samp$counts[,simplex]
    stats::dmultinom(counts,prob=p,log=TRUE)
}
ag2b <- function(ag,samp,c64=NULL){
    simplex <- c('Pb204','Pb206','Pb207','U238','UO2')
    tt <- hours(samp$time[,simplex])
    cc <- samp$counts[,simplex]
    dt <- samp$dwelltime[simplex]
    fn6 <- cc[,'Pb206'] ~ 1 + offset(ag['g']*tt[,'Pb206'])
    fn7 <- cc[,'Pb207'] ~ 1 + offset(ag['g']*tt[,'Pb207'])
    fnU <- cc[,'U238'] ~ 1 + offset(ag['g']*tt[,'U238'])
    fnUO <- cc[,'UO2'] ~ 1 + offset(ag['g']*tt[,'UO2'])
    b6 <- stats::glm(fn6,family=stats::poisson(link='log'))$coef
    b7 <- stats::glm(fn7,family=stats::poisson(link='log'))$coef
    bU <- stats::glm(fnU,family=stats::poisson(link='log'))$coef
    bUO <- stats::glm(fnUO,family=stats::poisson(link='log'))$coef
    if (is.null(c64)){
        b4 <- b6 - log(exp(ag['a4'])+1)
    } else {
        b4 <- b6 - log(c64*dt['Pb206']/dt['Pb204']) - log(exp(ag['a4'])+1)
    }
    out <- c(b4,b6,b7,bU,bUO)
    names(out) <- simplex
    out
}
b2p <- function(b,ag,tt){
    simplex <- c('Pb204','Pb206','Pb207','U238','UO2')
    nr <- nrow(tt)
    bm <- matrix(rep(b,nr),nrow=nr,byrow=TRUE)
    colnames(bm) <- simplex
    ebgt <- exp(bm + ag['g']*tt)
    out <- ebgt/sum(ebgt)
    colnames(out) <- simplex
    out
}

bias_correction <- function(samp,aLm,ag,c64=NULL){
    simplex <- c('Pb204','Pb206','Pb207','U238','UO2')
    tt <- hours(samp$time[,simplex])
    dt6 <- mean(tt[,'Pb206']-tt[,'U238'])
    dt7 <- mean(tt[,'Pb207']-tt[,'Pb206'])
    dtU <- mean(tt[,'UO2']-tt[,'U238'])
    L6c <- aLm$x['L6m'] - ag$x['g']*dt6
    L7c <- aLm$x['L7m'] - ag$x['g']*dt7
    LUc <- aLm$x['LUm'] - ag$x['g']*dtU
    if (is.null(c64)){
        L4c <- log(exp(ag$x['a4'])+1)
    } else {
        d4 <- samp$dwelltime['Pb204']
        d6 <- samp$dwelltime['Pb206']
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

get_gamma <- function(samp,ion='Pb206'){
    tt <- hours(samp$time[,ion])
    cc <- samp$counts[,ion]    
    fit <- stats::glm(cc ~ tt, family=stats::poisson(link='log'))
    fit$coef[2]
}

init_abO <- function(samp,oxide='UO2'){
    tt <- samp$time[,'U238']
    logUOU <- log(samp$cps[,oxide]) - log(samp$cps[,'U238'])
    fit <- lm(logUOU ~ tt)
    out <- fit$coef
    names(out) <- c('aO','bO')
    out
}
init_a6 <- function(samp,b6){
    tt <- samp$time[,'U238']
    logPbU <- log(samp$cps[,'Pb206']) - log(samp$cps[,'U238'])
    fit <- lm(logPbU ~ 1 + offset(b6*tt))
    fit$coef[1]
}

init_abg <- function(samp,A=0,B=1,oxide='UO2'){
    out <- rep(0,7)
    names(out) <- c('a4','a6','aO','b6','bO','g6','gO')
    out['g6'] <- get_gamma(samp,ion='Pb206')
    out['gO'] <- get_gamma(samp,ion=oxide)
    abO <- init_abO(samp,oxide=oxide)
    out['aO'] <- abO['aO']
    out['bO'] <- abO['bO']
    out['b6'] <- B*out['bO']
    out['a6'] <- init_a6(samp,b6=out['b6'])
    out['a4'] <- out['a6'] - A - B*out['aO']
    out
}
