logratios <- function(dat,c64=NULL){
    snames <- names(dat)
    ns <- length(snames)
    out <- list()
    for (sname in snames){
        print(sname)
        samp <- dat[[sname]]
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@2@"]]));##:ess-bp-end:##
        fm <- raw_count_ratios(samp=samp)
        afm <- avg_fm(fm)
        ag <- get_ag(samp=samp,c64=c64)
        out[[sname]] <- bias_correction(samp=samp,afm=afm,ag=ag)
    }
    out
}

raw_count_ratios <- function(samp){
    counts <- samp$counts[,c('206Pb','207Pb','238U','238U16O2')]
    lc6 <- log(counts[,'206Pb'])
    lc7 <- log(counts[,'207Pb'])
    lcU <- log(counts[,'238U'])
    lcUO <- log(counts[,'238U16O2'])
    nt <- length(lcU)
    f6m <- lc6 - lcU  # f6m = v in the paper
    f7m <- lc7 - lc6  # f7m = q in the paper
    fUm <- lcUO - lcU # fUm = w in the paper
    D <- (exp(f6m)+exp(f7m-f6m)+exp(fUm)+1)
    dDdf6m <- exp(f6m)-exp(f7m-f6m)
    dDdf7m <- exp(f7m-f6m)
    dDdfUm <- exp(fUm)
    d2Ddf6m2 <- exp(f6m)+exp(f7m-f6m)
    d2Ddf7m2 <- exp(f7m-f6m)
    d2DdfUm2 <- exp(fUm)
    d2Ddf6mdf7m <- exp(f7m-f6m)
    d2Ddf6mdfUm <- 0
    d2Ddf7mdfUm <- 0
    rs <- rowSums(counts)
    H11 <- (rs/D)*(dDdf6m*dDdf6m/D - d2Ddf6m2)
    H12 <- (rs/D)*(dDdf6m*dDdf7m/D - d2Ddf6mdf7m)
    H13 <- (rs/D)*(dDdf6m*dDdfUm/D - d2Ddf6mdfUm)
    H21 <- (rs/D)*(dDdf7m*dDdf6m/D - d2Ddf6mdf7m)
    H22 <- (rs/D)*(dDdf7m*dDdf7m/D - d2Ddf7m2)
    H23 <- (rs/D)*(dDdf7m*dDdfUm/D - d2Ddf7mdfUm)
    H31 <- (rs/D)*(dDdfUm*dDdf6m/D - d2Ddf6mdfUm)
    H32 <- (rs/D)*(dDdfUm*dDdf7m/D - d2Ddf7mdfUm)
    H33 <- (rs/D)*(dDdfUm*dDdfUm/D - d2Ddf7m2)
    out <- list()
    cnames <- c('f6m','f7m','fUm')
    for (i in 1:nt){
        out[[i]] <- list()
        out[[i]]$x <- c(f6m[i],f7m[i],fUm[i])
        H <- matrix(0,3,3)
        H[1,1] <- H11[i]; H[1,2] <- H12[i]; H[1,3] <- H13[i]
        H[2,1] <- H21[i]; H[2,2] <- H22[i]; H[2,3] <- H23[i]
        H[3,1] <- H31[i]; H[3,2] <- H32[i]; H[3,3] <- H33[i]
        out[[i]]$cov <- solve(-H)
        names(out[[i]]$x) <- cnames
        colnames(out[[i]]$cov) <- cnames
        rownames(out[[i]]$cov) <- cnames
    }
    out
}

avg_fm <- function(fm){
    nt <- length(fm$x)/2
    i1 <- 1:nt
    i2 <- (nt+1):(2*nt)
    f6m <- fm$x[i1]
    fUm <- fm$x[i2]
    E11 <- diag(fm$cov[i1,i1])
    E12 <- diag(fm$cov[i1,i2])
    E22 <- diag(fm$cov[i2,i2])
    adbc <- E11*E22-E12^2
    O11 <- E22/adbc
    O12 <- -E12/adbc
    O22 <- E11/adbc
    num6 <- sum(O22)*sum(f6m*O11+fUm*O12)-sum(O12)*sum(fUm*O22+f6m*O12)
    numU <- sum(O11)*sum(fUm*O22+f6m*O12)-sum(O12)*sum(f6m*O11+fUm*O12)
    den <- sum(O11)*sum(O22)-sum(O12)^2
    af6 <- num6/den
    afU <- numU/den
    out <- list()
    out$x <- c(af6,afU)
    Efm <- matrix(0,2,2)
    Efm[1,1] <- sum(O11)
    Efm[1,2] <- sum(O12)
    Efm[2,1] <- Efm[1,2]
    Efm[2,2] <- sum(O22)
    out$cov <- solve(Efm)
    labels <- c('f6m','fUm')
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

get_ag <- function(samp,c64=NULL){
    init <- init_ag(samp,c64=c64)
    np <- length(init) # number of parameters
    lower <- init - rep(1,np)
    upper <- init + rep(1,np)
    fit <- optim(init,misfit_ag,method='L-BFGS-B',
                 lower=lower,upper=upper,samp=samp,
                 c64=c64,control=list(fnscale=-1),hessian=TRUE)
    out <- list()
    out$x <- fit$par
    out$cov <- solve(-fit$hessian)
    out    
}
init_ag <- function(samp,c64=NULL){
    t4 <- hours(samp$time[,'204Pb'])
    tU <- hours(samp$time[,'238U'])
    c4 <- samp$counts[,'204Pb']
    cU <- samp$counts[,'238U']
    fitU <- glm(cU ~ tU, family=poisson(link='log'))
    bU <- fitU$coef[1]
    g <- fitU$coef[2]
    b4 <- glm(c4 ~ 1 + offset(g*t4), family=poisson(link='log'))$coef
    if (is.null(c64)){
        a4 <- log(exp(bU-b4)-1)
    } else {
        d4 <- samp$dwelltime['204Pb']
        d6 <- samp$dwelltime['206Pb']
        a4 <- log(exp(bU-b4)*(d4/d6)/c64 - 1)
    }
    out <- c(a4,g)
    names(out) <- c('a4','g')
    out
}
misfit_ag <- function(ag,samp,c64=NULL){
    simplex <- c('204Pb','206Pb','238U','238U16O2')
    b <- ag2b(ag,samp=samp,c64=c64)
    p <- b2p(b,ag['g'],tt=hours(samp$time[,simplex]))
    counts <- samp$counts[,simplex]
    dmultinom(counts,prob=p,log=TRUE)
}
ag2b <- function(ag,samp,c64=NULL){
    simplex <- c('204Pb','206Pb','238U','238U16O2')
    tt <- hours(samp$time[,simplex])
    cc <- samp$counts[,simplex]
    dt <- samp$dwelltime[simplex]
    fn6 <- cc[,'206Pb'] ~ 1 + offset(ag['g']*tt[,'206Pb'])
    fnU <- cc[,'238U'] ~ 1 + offset(ag['g']*tt[,'238U'])
    fnUO <- cc[,'238U16O2'] ~ 1 + offset(ag['g']*tt[,'238U16O2'])
    b6 <- glm(fn6,family=poisson(link='log'))$coef
    bU <- glm(fnU,family=poisson(link='log'))$coef
    bUO <- glm(fnUO,family=poisson(link='log'))$coef
    if (is.null(c64)){
        b4 <- b6 - log(exp(ag['a4'])+1)
    } else {
        b4 <- b6 - log(c64*dt['206Pb']/dt['204Pb']) - log(exp(ag['a4'])+1)
    }
    out <- c(b4,b6,bU,bUO)
    names(out) <- simplex
    out
}
b2p <- function(b,ag,tt){
    simplex <- c('204Pb','206Pb','238U','238U16O2')
    nr <- nrow(tt)
    bm <- matrix(rep(b,nr),nrow=nr,byrow=TRUE)
    colnames(bm) <- simplex
    ebgt <- exp(bm + ag['g']*tt - bm[,'238U'])
    out <- ebgt/sum(ebgt)
    colnames(out) <- simplex
    out
}

bias_correction <- function(samp,afm,ag,c64=NULL){
    simplex <- c('204Pb','206Pb','238U','238U16O2')
    tt <- hours(samp$time[,simplex])
    dt <- colMeans(tt[,c(1,2,4)]-tt[,3])
    names(dt) <- c('206Pb','238U16O2')
    f6c <- afm$x['f6m'] - ag$x['g']*dt['206Pb']
    fUc <- afm$x['fUm'] - ag$x['g']*dt['238U16O2']
    if (is.null(c64)){
        f4c <- f6c - log(exp(ag$x['a4'])+1)
    } else {
        d4 <- samp$dwelltime['204Pb']
        d6 <- samp$dwelltime['206Pb']
        f4c <- f6c - log(exp(ag$x['a4'])+1) - log(c64*d6/d4)
    }
    out <- list()
    out$x <- c(f4c,f6c,fUc)
    # joint covariance matrix of afm and ag:
    E <- matrix(0,4,4)
    E[1:2,1:2] <- afm$cov
    E[3:4,3:4] <- ag$cov
    J <- matrix(0,3,4)
    J[1,1] <- 1 # df4c/df6m
    J[2,1] <- 1 # df6c/df6m
    J[3,2] <- 1 # dfUc/dfUm
    J[1,3] <- -exp(ag$x['a4'])/(exp(ag$x['a4'])+1) # df4c/da4
    J[1,4] <- -dt['206Pb'] # df4cdg
    J[2,4] <- -dt['206Pb'] # df6cdg
    J[3,4] <- -dt['238U16O2'] # dfUcdg
    out$cov <- J %*% E %*% t(J)
    labels <- c('f4c','f6c','fUc')
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

lr2XY <- function(lr,dat,c64=NULL){
    snames <- names(dat)
    ns <- length(snames)
    out <- matrix(0,ns,5)
    colnames(out) <- c('X','sX','Y','sY','rho')
    rownames(out) <- snames
    J <- matrix(0,2,3)
    rownames(J) <- c('X','Y')
    colnames(J) <- c('f4c','f6c','fUc')
    J['X','fUc'] <- 1
    for (sname in snames){
        d4 <- dat[[sname]]$dwelltime['204Pb']
        d6 <- dat[[sname]]$dwelltime['206Pb']
        dU <- dat[[sname]]$dwelltime['238U']
        dUO <- dat[[sname]]$dwelltime['238U16O2']
        f4c <- lr[[sname]]$x['f4c']
        f6c <- lr[[sname]]$x['f6c']
        fUc <- lr[[sname]]$x['fUc']
        X <- fUc + log(dU/dUO)
        if (is.null(c64)){
            Y <- f6c + log(dU/d6)
            J['Y','f6c'] <- 1
        } else {
            Y <- log(dU) + log(exp(f6c)/d6 - exp(f4c)*c64/d4)
            den <- exp(f6c)*d4 - exp(f4c)*c64*d6
            J['Y','f4c'] <- -exp(f4c)*c64*d6/den
            J['Y','f6c'] <- exp(f6c)*d4/den
        }
        E <- J %*% lr[[sname]]$cov %*% t(J)
        out[sname,'X'] <- X
        out[sname,'Y'] <- Y
        out[sname,'sX'] <- sqrt(E['X','X'])
        out[sname,'sY'] <- sqrt(E['Y','Y'])
        out[sname,'rho'] <- cov2cor(E)['X','Y']
    }
    out
}

# build a calibration curve
calibration <- function(XY,plot=TRUE,disp=TRUE){
    xlab <- quote(''^238*'U'^16*'O'[2]*'/'^238*'U')
    ylab <- quote(''^206*'Pb*/'^238*'U')
    fit <- IsoplotR::isochron(XY,plot=plot,
                              xlab=xlab,ylab=ylab)
    out <- list()
    out$x <- c(fit$a['a'],fit$b['b'])
    out$cov <- matrix(0,2,2)
    if (disp) d <- fit$mswd
    else d <- 1
    out$cov[1,1] <- d*fit$a['s[a]']^2
    out$cov[2,2] <- d*fit$b['s[b]']^2
    out$cov[1,2] <- d*fit$cov.ab
    out$cov[2,1] <- out$cov[1,2]
    labels <- c('A','B')
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

# apply a calibration to a unknown samples
calibrate <- function(XY,fit,PbUstand=NULL,tst=c(337.13,0.18)){
    ns <- nrow(XY)
    if (is.null(PbUstand))
        Pb6U8stand <- IsoplotR:::age_to_Pb206U238_ratio(tt=tst[1],st=tst[2])
    Ytstand <- log(Pb6U8stand[1])
    out <- list()
    out$x <- XY[,'Y'] + Ytstand - fit$x['A'] - fit$x['B']*XY[,'X']
    sYtstand <- Pb6U8stand[2]/Pb6U8stand[1]
    E <- matrix(0,2*ns+3,2*ns+3)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    E[i1,i1] <- diag(XY[,'sX']^2)
    E[i2,i2] <- diag(XY[,'sY']^2)
    E[i1,i2] <- diag(XY[,'rho']*XY[,'sX']*XY[,'sY'])
    E[i2,i1] <- E[i1,i2]
    E[(2*ns+1),(2*ns+1)] <- sYtstand^2
    E[(2*ns)+(2:3),(2*ns)+(2:3)] <- fit$cov
    J <- matrix(0,nrow=ns,ncol=2*ns+3)
    rownames(J) <- rownames(XY)
    J[i1,i1] <- diag(-fit$x['B'],ns,ns)
    J[i1,i2] <- diag(1,ns,ns)
    J[,2*ns+1] <- 1
    J[,2*ns+2] <- -1
    J[,2*ns+3] <- -XY[,'X']
    out$cov <- J %*% E %*% t(J)
    out
}

ages <- function(logPbU){
    out <- list()
    Pb6U8 <- exp(logPbU$x)
    JR <- diag(Pb6U8)
    E <- JR %*% logPbU$cov %*% t(JR)
    l38 <- IsoplotR::settings('lambda','U238')[1]
    out$x <- log(1+Pb6U8)/l38
    Jt <- diag(1/(l38*(1+Pb6U8)))
    out$cov <- Jt %*% E %*% t(Jt)
    labels <- names(logPbU$x)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

predict_cps <- function(samp,c64=NULL){
    simplex <- c('204Pb','206Pb','238U','238U16O2')
    fm <- raw_count_ratios(samp=samp)
    afm <- avg_fm(fm)
    eafm <- exp(afm$x)
    ag <- get_ag(samp=samp)
    d4 <- samp$dwelltime['204Pb']
    d6 <- samp$dwelltime['206Pb']
    if (is.null(c64)){
        f4m <- afm$x['f6m'] - log(exp(ag$x['a4'])+1)
    } else {
        f4m <- afm$x['f6m'] - log(exp(ag$x['a4'])+1) - log(c64*d6/d4)
    }
    den <- sum(eafm) + 1
    p4 <- exp(f4m)/den
    p6 <- eafm['f6m']/den
    pU <- 1/den
    pUO <- eafm['fUm']/den
    p <- c(p4,p6,pU,pUO)
    rs <- rowSums(samp$counts[,simplex])
    dt <- samp$dwelltime[simplex]
    counts <- matrix(rs,ncol=1) %*% matrix(p,nrow=1)
    cps <- sweep(counts,MARGIN=2,FUN='/',dt)
    colnames(cps) <- simplex
    cps
}

hours <- function(tt){
    tt/3600
}
