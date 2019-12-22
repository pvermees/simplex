#' @title apply a calibration to sample data
#' @description fit a straight line to log[Pb/U] vs. log[UOx/U] sample
#'     data
#' @details Regresses a line through
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U and
#'     Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U using the same
#'     slope as previously determined on an age standard
#' @param dat a dataset of class \code{simplex}
#' @param fit the output of \code{calibration}
#' @param syserr propagate the systematic error associated with the
#'     standard (age uncertainty and calibration fit)?
#' @return an \code{IsoplotR} object of class \code{UPb}
#'     (\code{format=5})
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' fit <- calibration(stand,oxide='UO2')
#' unk <- unknowns(dat=Cameca,prefix='Plesovice',invert=TRUE)
#' cal <- calibrate(unk,fit)
#' @export
calibrate <- function(dat,fit,syserr=FALSE){
    # 1. calibrate, calculate (co)variances and partial derivatives
    snames <- names(dat)
    ns <- length(snames)
    X <- matrix(0,ns,8)
    rownames(X) <- snames
    for (sname in snames){
        calspot <- calibrate_spot(dat[[sname]],fit)
        X[sname,] <- calspot
    }
    colnames(X) <- names(calspot)
    # 2. Collate calibrated data into one big covariance structure
    out <- list()
    class(out) <- 'UPb'
    out$format <- 9
    out$names <- rownames(X)
    out$logratios <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
    ns <- length(out$names)
    nr <- length(out$logratios)
    out$x <- rep(0,ns*nr)
    covmat <- matrix(0,ns*nr,ns*nr)
    J <- matrix(0,nrow=ns*nr,ncol=3)
    colnames(J) <- c('A','B','U8Pb6s')
    U8Pb6s <- -log(fit$PbU[1])
    varU8Pb6s <- fit$PbU[2]/fit$PbU[1]
    for (i in 1:ns){
        j <- (i-1)*nr+(1:nr)
        out$x[j[1]] <- U8Pb6s + fit$AB['A'] - X[i,'Ax']
        out$x[j[2:3]] <- X[i,c('Pb76','Pb46')]
        J[j[1],'A'] <- 1
        J[j[1],'B'] <- -X[i,'dAxdBs']
        J[j[1],'U8Pb6s'] <- 1
        covmat[j[1],j[1]] <- X[i,'varAx']
        covmat[j[2],j[2]] <- X[i,'varPb76']
        covmat[j[3],j[3]] <- X[i,'varPb46']
        covmat[j[2],j[3]] <- X[i,'covPb46Pb76']
        covmat[j[3],j[2]] <- covmat[j[2],j[3]]
    }
    E <- matrix(0,3,3)
    E[1:2,1:2] <- fit$cov
    E[3,3] <- varU8Pb6s
    if (syserr) out$cov <- covmat + J %*% E %*% t(J)
    else out$cov <- covmat
    out
}

calibrate_spot <- function(spot,fit){
    out <- rep(0,8)
    names(out) <- c('Ax','Pb76','Pb46',
                    'varAx','varPb76','varPb46','covPb46Pb76',
                    'dAxdBs')
    p <- pars(spot,oxide=fit$oxide)
    g <- get_gamma(B=fit$AB['B'],p=p)
    init <- log(sum(p$c6))-log(sum(p$cU))
    Ax <- stats::optimise(misfit_A,interval=c(-10,10),
                          spot=spot,fit=fit)$minimum
    out['Ax'] <- Ax
    out['varAx'] <- solve(stats::optimHess(Ax,misfit_A,spot=spot,fit=fit))
    out['dAxdBs'] <- dAxdBs(p,g,Ax=Ax,Bs=fit$AB['B'])
    out['Pb76'] <- getPbLogRatio(g=g$Pb,tn=p$t7,td=p$t6,
                                 nn=p$n7,nd=p$n6,dn=p$d7,dd=p$d6)
    out['Pb46'] <- getPbLogRatio(g=g$Pb,tn=p$t4,td=p$t6,
                                 nn=p$n4,nd=p$n6,dn=p$d4,dd=p$d6)
    HPb <- matrix(0,2,2)
    sumn <- sum(p$n4+p$n6+p$n7)
    HPb[1,1] <- -sum(p$n4)*sum(p$n6+p$n7)/sumn
    HPb[1,2] <- sum(p$n4)*sum(p$n7)/sumn
    HPb[2,2] <- -sum(p$n7)*sum(p$n6+p$n4)/sumn
    HPb[2,1] <- HPb[1,2]/sumn
    covmat <- solve(-HPb)
    out['varPb46'] <- covmat[1,1]
    out['varPb76'] <- covmat[2,2]
    out['covPb46Pb76'] <- covmat[1,2]
    out
}

getPbLogRatio <- function(g,tn,td,nn,nd,dn,dd){
    init <- log(sum(nn/dn)) - log(sum(nd/dd)) + g*mean(td-tn)
    interval <- c(0.9*init,1.1*init)
    optimise(LL_Pb,interval=interval,g=g,tn=tn,td=td,
             nn=nn,nd=nd,dn=dn,dd=dd,maximum=TRUE)$maximum
}

LL_Pb <- function(b,g,tn,td,nn,nd,dn,dd){
    bn<- b + g*(tn-td) + log(dn)-log(dd)
    LL <- nn*bn - (nn+nd)*log(1+exp(bn))
    sum(LL)
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
