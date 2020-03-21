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
        if (syserr) J[i,ns+1] <- calspot['dAdB']
    }
    EAx <- J %*% EAxB %*% t(J)
    # 2. get PD ratios
    out <- list()
    out$snames <- snames
    out$num <- fit$daughter
    out$den <- fit$parent
    dp <- paste0(fit$daughter,fit$parent)
    out$x <- log(fit$DP[dp]) + Ax - fit$AB['A']
    EAxAslDPs <- matrix(0,ns+2,ns+2)
    EAxAslDPs[1:ns,1:ns] <- EAx
    EAxAslDPs[(ns+1),(ns+1)] <- fit$AB.cov['A','A']
    EAxAslDPs[(ns+2),(ns+2)] <- fit$DP.cov[dp,dp]
    J <- matrix(0,ns,ns+2)
    J[(1:ns),(1:ns)] <- -diag(ns)        # dPDx_dAx
    if (syserr) {
        J[(1:ns),(ns+1)] <- 1            # dPDx_dAs
        J[(1:ns),(ns+2)] <- 1/fit$DP[dp] # dPDx_dDPs
    }
    out$cov <- J %*% EAxAslDPs %*% t(J)
    out
}

calibrate_spot <- function(B,p,b0g){
    out <- rep(0,3)
    names(out) <- c('A','varA','dAdB')
    x <- log(sum(p[[p$oxide]]$c)) - log(sum(p[[p$parent]]$c))
    y <- log(sum(p[[p$daughter]]$c)) - log(sum(p[[p$parent]]$c))
    init <- y - B * x
    out['A'] <- stats::optimise(misfit_A,interval=init+c(-5,5),
                                 B=B,p=p,b0g=b0g)$minimum
    H <- stats::optimHess(par=out['A'],fn=misfit_A,B=B,p=p,b0g=b0g)
    out['varA'] <- solve(H)
    out['dAdB'] <- misfit_A(A=out['A'],B=B,p=p,b0g=b0g,deriv=TRUE)
    out
}

mergecal <- function(...){
    cals <- list(...)
    nc <- length(cals)
    out <- list()
    cal <- cals[[1]]
    out$num <-cal$num
    out$den <-cal$den
    out$snames <- cal$snames
    out$x <- cal$x
    out$cov <- cal$cov
    if (nc<2) return(out)
    for (i in 2:nc){
        cal <- cals[[i]]
        if (!all(out$snames%in%cal$snames))
            stop('Mismatched sample names.')
        nold <- length(out$x)
        nnew <- length(cal$x)
        out$num <- c(out$num,cal$num)
        out$den <- c(out$den,cal$den)
        out$x <- c(out$x,cal$x)
        covmat <- diag(length(out$x))
        covmat[1:nold,1:nold] <- out$cov
        covmat[(nold+1):(nold+nnew),(nold+1):(nold+nnew)] <- cal$cov
        out$cov <- covmat
    }
    out
}
