# get geometric mean Pb207/Pb206 ratio to estimate
# the standard age if not supplied by the user
get_Pb76 <- function(dat){
    snames <- names(dat)
    ns <- length(snames)
    lPb76 <- rep(0,ns)
    for (i in 1:ns){
        p <- pars(dat[[i]])
        lPb76[i] <- log(sum(p$c7)/sum(p$c6))
    }
    lPb76 <- mean(lPb76)
    slPb76 <- stats::sd(lPb76)/sqrt(ns)
    Pb76 <- exp(lPb76)
    sPb76 <- Pb76*slPb76
    c(Pb76,sPb76)
}

get_gamma <- function(B,p){
    fitU <- optim(c(log(p$cU[1]),0),LL_gamma,nn=p$nU,dd=p$dU,tt=p$tU)
    fitO <- optim(c(log(p$cO[1]),0),LL_gamma,nn=p$nO,dd=p$dO,tt=p$tO)
    out <- list()
    out$U <- fitU$par[2]
    out$O <- fitO$par[2]
    out$Pb <- B*out$O + (1-B)*out$U
    out
}
LL_gamma <- function(pars,nn,dd,tt){
    intercept <- pars[1]
    slope <- pars[2]
    LL <- nn*(intercept + slope*tt + log(dd)) - exp(intercept + slope*tt)*dd
    -sum(LL)
}

getPbLogRatio <- function(p,g,num,c64=1){
    ax <- get_alpha(p=p,g=g,den=num,c64=c64)
    a6 <- get_alpha(p=p,g=g,den=6,c64=c64)
    optimise(LL_Pb,interval=c(-20,20),p=p,g=g,
             anum=ax,aden=a6,num=num,maximum=TRUE)$maximum
}
LL_Pb <- function(b,p,g,anum,aden,num){
    tx <- p[[paste0('t',num)]]
    dx <- p[[paste0('d',num)]]
    nx <- p[[paste0('n',num)]]
    log_nx6i <- b + anum - aden - log(p$d6/dx) - g*(p$t6-tx)
    LL <- nx*log_nx6i - (nx+p$n6)*log(1+exp(log_nx6i))
    sum(LL)
}

get_alpha <- function(p,g,den,c64=1){
    tx <- p[[paste0('t',den)]]
    maxb <- - max(g*tx) - 1e-5 - log(c64)
    b <- optimise(LL_bb,p=p,g=g,den=den,lower=-10,upper=maxb,
                  maximum=TRUE)$maximum
    log(1-exp(b+g*tx))
}
LL_bb <- function(b,p,g,den){
    tx <- p[[paste0('t',den)]]
    dx <- p[[paste0('d',den)]]
    nx <- p[[paste0('n',den)]]
    bdx <- p[[paste0('bd',den)]]
    bnx <- p[[paste0('bn',den)]]
    log_nbx <- b + g*tx + log(bdx/dx)
    LL <- bnx*log_nbx - (bnx+nx)*log(1+exp(log_nbx))
    sum(LL)
}

logratios2ratios <- function(lr){
    out <- list()
    out$x <- exp(lr$x)
    J <- diag(out$x)
    out$cov <- J %*% lr$cov %*% t(J)
    out
}
