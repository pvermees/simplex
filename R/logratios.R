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
    # slope(U), intercept(U), slope(O), intercept(O), intercept(Pb) :
    init <- c(log(p$cU[1]),0,log(p$cO[1]),0,log(p$c6[1]))
    fit <- optim(par=init,fn=LL_gamma,B=B,p=p)$par
    out <- list()
    out$U <- fit[2]
    out$O <- fit[4]
    out$Pb <- B*fit[4] + (1-B)*fit[2]
    out
}
LL_gamma <- function(pars,B,p){
    log_nU <- pars[1] + pars[2]*p$tU
    LL_U <- p$nU*(log_nU + log(p$dU)) - exp(log_nU)*p$dU
    log_nO <- pars[3] + pars[4]*p$tO
    LL_O <- p$nO*(log_nO + log(p$dO)) - exp(log_nO)*p$dO
    pars[6] <- B*pars[4] + (1-B)*pars[2]
    log_n6 <- pars[5] + pars[6]*p$t6
    LL_6 <- p$n6*(log_n6 + log(p$d6)) - exp(log_n6)*p$d6
    -sum(LL_U + LL_O + LL_6)
}

getPbLogRatio <- function(p,g,num,c64=NULL){
    ax <- get_alpha(p=p,g=g,den=num)
    a6 <- get_alpha(p=p,g=g,den=6)
    if (is.null(c64)) interval <- c(-20,20)
    else interval <- -log(c64)-c(20,1e-10)
    optimise(LL_beta,p=p,g=g,ax=ax,a6=a6,num=num,
             interval=interval,maximum=TRUE)$maximum
}
LL_beta <- function(b,p,g,ax,a6,num){
    tx <- p[[paste0('t',num)]]
    dx <- p[[paste0('d',num)]]
    nx <- p[[paste0('n',num)]]
    bi <- b + g*(tx-p$t6) + ax - a6
    log_nx6i <- bi + log(dx/p$d6)
    LL <- nx*log_nx6i -(nx+p$n6)*log(1+exp(log_nx6i))
    sum(LL)
}

get_alpha <- function(p,g,den){
    tx <- p[[paste0('t',den)]]
    interval <- c(-10,0)-(max(g*tx)+1e-5)
    b <- optimise(LL_balpha,p=p,g=g,den=den,
                  interval=interval,maximum=TRUE)$maximum
    log(1-exp(b+g*tx))
}
LL_balpha <- function(b,p,g,den){
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
