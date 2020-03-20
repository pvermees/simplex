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
    # 1. get sample intercepts
    snames <- names(dat)
    ns <- length(snames)
    Ax <- rep(0,ns)
    EAxB <- matrix(0,ns+1,ns+1)
    EAxB[(ns+1),(ns+1)] <- fit$AB.cov['B','B']
    J <- matrix(0,ns,ns+1)
    J[(1:ns),(1:ns)] <- diag(ns)
    for (i in 1:ns){
        spot <- dat[[i]]
        p <- pars(spot=spot,oxide=fit$oxide)
        b0g <- get_b0g(spot=spot,oxide=fit$oxide)
        calspot <- calibrate_spot(B=fit$AB['B'],p=p,b0g=b0g)
        Ax[i] <- calspot['A']
        EAxB[i,i] <- calspot['varA']
        J[i,ns+1] <- calspot['dAdB']
    }
    EAx <- J %*% EAxB %*% t(J)
    # 2. get PD ratios
    out <- list()
    out$parent <- fit$parent
    out$daughter <- fit$daughter
    if (fit$parent=='U238') dp <- 'Pb206U238'
    else if (fit$parent=='Th232') dp <- 'Pb208Th232'
    else stop('Illegal parent used in calibration.')
    out$PD <- -(log(fit$DP[dp]) + fit$AB['A'] - Ax)
    EAxAsDPs <- matrix(0,ns+2,ns+2)
    EAxAsDPs[1:ns,1:ns] <- EAx
    EAxAsDPs[(ns+1),(ns+1)] <- fit$AB.cov['A','A']
    EAxAsDPs[(ns+2),(ns+2)] <- fit$DP.cov[dp,dp]
    J <- matrix(0,ns,ns+2)
    J[(1:ns),(1:ns)] <- diag(ns)      # dPDx_dAx
    J[(1:ns),(ns+1)] <- -1            # dPDx_dAs
    J[(1:ns),(ns+2)] <- -1/fit$DP[dp] # dPDx_dDPs
    out$PD.cov <- J %*% EAxAsDPs %*% t(J)
    names(out$PD) <- snames
    rownames(out$PD.cov) <- snames
    colnames(out$PD.cov) <- snames
    out
}

calibrate_spot <- function(B,p,b0g){
    out <- rep(0,3)
    names(out) <- c('A','varA','dAdB')
    x <- log(sum(p$O$c)) - log(sum(p$U238$c))
    y <- log(sum(p$Pb206$c)) - log(sum(p$U238$c))
    init <- y - B * x
    out['A'] <- stats::optimise(misfit_A,interval=init+c(-5,5),
                                 B=B,p=p,b0g=b0g)$minimum
    H <- stats::optimHess(par=out['A'],fn=misfit_A,B=B,p=p,b0g=b0g)
    out['varA'] <- solve(H)
    out['dAdB'] <- misfit_A(A=out['A'],B=B,p=p,b0g=b0g,deriv=TRUE)
    out
}
