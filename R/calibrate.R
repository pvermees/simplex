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
    snames <- names(dat)
    ns <- length(snames)
    B0G <- get_B0G(dat=dat,parent=fit$parent,oxide=fit$oxide)
    misfit <- function(par){
        out <- 0
        for (i in 1:ns){
            spot <- dat[[i]]
            p <- pars(spot=spot,oxide=fit$oxide)
            As <- fit$AB['A'] + par[i]
            out <- out + misfit_A(As,B=fit$AB['B'],p=p,b0g=B0G[[i]])
        }
        out
    }
    init <- rep(0,ns)
    for (i in 1:ns){
        spot <- dat[[i]]
        p <- pars(spot=spot,oxide=fit$oxide)
        x <- log(sum(p[[p$oxide]]$c)) - log(sum(p[[p$parent]]$c))
        y <- log(sum(p[[p$daughter]]$c)) - log(sum(p[[p$parent]]$c))
        init[i] <- y - fit$AB['B'] * x - fit$AB['A']
    }
    out <- list()
    fit <- stats::optim(init,misfit,method='L-BFGS-B',
                        lower=init-1,upper=init+1,hessian=TRUE)
    out$x <- fit$par
    out$cov <- solve(fit$hessian)
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
