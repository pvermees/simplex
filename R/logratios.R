logratios <- function(dat,fit){
    X <- calibrate(dat=dat,fit=fit)
    snames <- names(dat)
    ns <- length(snames)
    ratios <- c('Pb206U238','Pb207Pb206','Pb204Pb206')
    nr <- length(ratios)
    out <- list()
    out$snames <- snames
    out$x <- rep(0,ns*nr)
    covmat <- matrix(0,ns*nr,ns*nr)
    J <- matrix(0,nrow=ns*nr,ncol=3)
    colnames(J) <- c('A','B','Pb6U8s')
    for (i in 1:ns){
        j <- (i-1)*nr+(1:nr)
        out$x[j] <- X[i,c('Pb6U8','Pb76','Pb46')]
        J[j[1],'A'] <- -1
        J[j[1],'B'] <- X[i,'dPb6U8dBs']
        J[j[1],'Pb6U8s'] <- 1
        covmat[j[1],j[1]] <- X[i,'varPb6U8']
        covmat[j[2],j[2]] <- X[i,'varPb76']
        covmat[j[3],j[3]] <- X[i,'varPb46']
        covmat[j[2],j[3]] <- X[i,'covPb46Pb76']
        covmat[j[3],j[2]] <- covmat[j[2],j[3]]
    }
    E <- matrix(0,3,3)
    E[1:2,1:2] <- fit$cov
    out$cov <- covmat + J %*% E %*% t(J)
    out
}

calibrate <- function(dat,fit){
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
    Ax <- optimise(misfit_A,interval=c(-10,10),spot=spot,fit=fit)$minimum
    As <- fit$AB['A']
    Bs <- fit$AB['B']
    out['Pb6U8'] <- log(fit$PbU[1]) + Ax - As
    out['varPb6U8'] <- solve(optimHess(Ax,misfit_A,spot=spot,fit=fit))
    out['dPb6U8dBs'] <- dAxdBs(p,g,Ax=Ax,Bs=Bs)
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

misfit_A <- function(A,spot,fit){
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
    LL <- sum(counts*log(theta))
    -LL
}

# implicit differentiation of dLLdAx to get dAxdBs
dAxdBs <- function(p,g,Ax,Bs){
    aO <- log(p$cO/p$cU) + g$O*(p$tU-p$tO)
    n6U <- exp( Ax + Bs*aO + g$Pb*(p$t6-p$tU) )*p$d6/p$dU
    dn6UdAx <- n6U
    dn6UdBs <- n6U*aO
    nt <- length(p$tU)
    theta6 <- n6U/(n6U+1)
    thetaU <- 1/(n6U+1)
    dtheta6dAx <- n6U/(n6U+1)^2
    dthetaUdAx <- -n6U/(n6U+1)^2
    dtheta6dBs <- aO*dtheta6dAx
    dthetaUdBs <- aO*dthetaUdAx
    d2theta6dAx2 <- (n6U-1)/(n6U+1)
    d2thetaUdAx2 <- (1-n6U)/(n6U+1)
    d2theta6dAxdBs <- aO*(n6U-1)/(n6U+1)
    d2thetaUdAxdBs <- aO*(1-n6U)/(n6U+1)
    d2LLdAxdBs <- sum( (p$n6/theta6)*d2theta6dAxdBs +
                       (p$nU/thetaU)*d2thetaUdAxdBs -
                       (p$n6/theta6^2)*dtheta6dAx*dtheta6dBs -
                       (p$nU/thetaU^2)*dthetaUdAx*dthetaUdBs )
    d2LLdAx2 <- sum( (p$n6/theta6)*d2theta6dAx2 +
                     (p$nU/thetaU)*d2thetaUdAx2 -
                     (p$n6/theta6^2)*dtheta6dAx^2 -
                     (p$nU/thetaU^2)*dthetaUdAx^2 )
    -d2LLdAxdBs/d2LLdAx2
}
