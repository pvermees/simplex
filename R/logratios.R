logratios <- function(dat,fit){
    snames <- names(dat)
    ns <- length(snames)
    out <- matrix(0,ns,8)
    dPb6U8dAs <- -1
    rownames(out) <- snames
    for (sname in snames){
        calspot <- calibrate_spot(dat[[sname]],fit)
        out[sname,] <- calspot
    }
    colnames(out) <- names(calspot)
    out
}

calibrate_spot <- function(spot,fit){
    out <- rep(0,8)
    names(out) <- c('Pb6U8','Pb76','Pb46',
                    'varPb6U8','varPb76','varPb46','covPb46Pb76',
                    'dPb6U8dBs')
    p <- pars(spot,oxide=fit$oxide)
    g <- get_gamma(B=fit$AB['B'],p=p)
    init <- log(sum(p$c6))-log(sum(p$cU))
    A <- optimise(LL_A,interval=c(-10,10),spot=spot,fit=fit)$minimum
    vA <- solve(optimHess(A,LL_A,spot=spot,fit=fit))
    out['Pb6U8'] <- log(fit$PbU[1]) + A - fit$AB['A']
    out['varPb6U8'] <- vA
    out['dPb6U8dBs'] <- 0 # TODO
    out['Pb46'] <- log(0.5+sum(p$c4))-log(0.5+sum(p$c6))+g$Pb*mean(p$t6-p$t4)
    out['Pb76'] <- log(sum(p$c7))-log(sum(p$c6))+g$Pb*mean(p$t6-p$t7)
    HPb <- matrix(0,2,2)
    nsum <- sum(p$n4+p$n6+p$n7)
    HPb[1,1] <- -sum(p$n4)*sum(p$n6+p$n7)/nsum
    HPb[1,2] <- sum(p$n4)*sum(p$n7)/nsum
    HPb[2,2] <- -sum(p$n7)*sum(p$n6+p$n4)/nsum
    HPb[2,1] <- HPb[1,2]/nsum
    covmat <- solve(-HPb)
    out['varPb46'] <- covmat[1,1]
    out['varPb76'] <- covmat[2,2]
    out['covPb46Pb76'] <- covmat[1,2]
    out
}

LL_A <- function(A,spot,fit){
    B <- fit$AB['B']
    AB <- c(A,B)
    p <- pars(spot,oxide=fit$oxide)
    g <- get_gamma(B=B,p=p)
    aO <- log(p$cO/p$cU) + g$O*(p$tU-p$tO)
    a6 <- A + B*aO
    # 1. get count ratios
    n6U <- exp(a6 + g$Pb*(p$t6-p$tU))*p$d6/p$dU
    # 2. get proportions
    nt <- length(p$tU)
    den <- n6U + 1
    theta <- cbind(n6U,1)/matrix(rep(den,2),ncol=2)
    colnames(theta) <- c('Pb206','U238')
    counts <- spot$counts[,c('Pb206','U238')]
    # 3. compute multinomial log-likelihood
    out <- stats::dmultinom(x=counts,prob=theta,log=TRUE)
    -out
}
