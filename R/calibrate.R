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
    init <- rep(0,ns)
    for (i in 1:ns){
        spot <- dat[[i]]
        p <- pars(spot=spot,oxide=fit$oxide)
        x <- log(sum(p[[p$oxide]]$c)) - log(sum(p[[p$parent]]$c))
        y <- log(sum(p[[p$daughter]]$c)) - log(sum(p[[p$parent]]$c))
        init[i] <- y - fit$AB['B'] * x - fit$AB['A']
    }
    B0G <- get_B0G(dat=dat,parent=fit$parent,oxide=fit$oxide)
    # par is the vertical difference between the sample and the standard curve
    misfit_spot <- function(par,spot,fit,b0g,deriv=FALSE){
        p <- pars(spot=spot,parent=fit$parent,
                  daughter=fit$daughter,oxide=fit$oxide)
        XY <- getCalXY(p,b0g)
        Yp <- fit$AB['A'] + fit$AB['B']*XY$X + par
        bDPpc <- Yp + A2Corr(p=p,b0g=b0g,num=p$daughter,den=p$parent)
        bn <- bDPpc + log(p[[p$daughter]]$d) - log(p[[p$parent]]$d)
        nnum <- p[[p$daughter]]$n
        nden <- p[[p$parent]]$n
        if (deriv){
            dLL_dbn <- nnum - (nnum+nden)*exp(bn)/(1+exp(bn))
            dbn_dpar <- 1
            dbn_dA <- 1
            dbn_dB <- XY$X
            dLL_dpar <- sum(dLL_dbn*dbn_dpar)
            dLL_dA <- sum(dLL_dbn*dbn_dA)
            dLL_dB <- sum(dLL_dbn*dbn_dB)
            dpar_dA <- -dLL_dA/dLL_dpar
            dpar_dB <- -dLL_dB/dLL_dpar
            out <- -c(dpar_dA,dpar_dB)
        } else {
            out <- -sum(LLbinom(bn=bn,nnum=nnum,nden=nden))
        }
        out
    }
    misfit <- function(par,dat,fit,B0G){
        out <- 0
        for (i in 1:ns){
            out <- out + misfit_spot(par=par[i],spot=dat[[i]],
                                     fit=fit,b0g=B0G[[i]])
        }
        out
    }
    jacobian <- function(par,dat,fit,B0G){
        ns <- length(dat)
        out <- matrix(0,ns,2)
        for (i in 1:ns){
            out[i,] <- misfit_spot(par=par[i],spot=dat[[i]],
                                   fit=fit,b0g=B0G[[i]],deriv=TRUE)
        }
        out        
    }
    dAfit <- stats::optim(init,misfit,method='BFGS',hessian=TRUE,
                          dat=dat,fit=fit,B0G=B0G)
    dp <- paste0(fit$daughter,fit$parent)
    out <- fit
    class(out) <- "calibrated"
    out$num <- fit$daughter
    out$den <- fit$parent
    out$snames <- snames
    out$x <- log(fit$DP[dp]) + dAfit$par
    out$cov <- solve(dAfit$hessian)
    if (syserr){
        J <- jacobian(par=dAfit$par,dat=dat,fit=fit,B0G=B0G)
        E <- fit$AB.cov
        out$cov <- out$cov + J %*% E %*% t(J)
    }
    out
}

#' @title merge calibrated data
#' @description groups calibrated U-Pb, Th-Pb and Pb-Pb ratios into
#'     one object for further processing and plotting
#' @param ... any number of objects of class \code{calibrated}
#' @return a new object of class \code{calibrated}
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' calU <- calibration(stand=stand,oxide='UO2',parent='U238',
#'                     daughter='Pb206',cD4=18.7)
#' samp <- unknowns(dat=Cameca,prefix='Qinghu')
#' calUsamp <- calibrate(dat=samp,fit=calU)
#' calPbsamp <- getPbLogRatios(samp)
#' calsamp <- mergecal(calUsamp,calPbsamp)
#' tab <- data2table(calsamp)
#' @export
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
    class(out) <- "calibrated"
    out
}
