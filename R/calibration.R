#' @title create a calibration curve
#' @description fit a straight line to log[Pb/U] vs. log[UOx/U] data
#' @details Regresses a line through
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U and
#'     Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U
#' @param stand a dataset of class \code{simplex}
#' @param oxide either \code{'UO'} or \code{'UO2'}
#' @param omit indices of measurements to be omitted from the
#'     calibration
#' @return a list with the following items:
#' 
#' \code{AB} a vector with the intercept (\code{'A'}) and slope
#' (\code{'B'}) of the power law
#'
#' \code{cov} the covariance matrix of \code{AB}
#'
#' @examples
#' data(Cameca,package="simplex")
#' stand <- standards(dat=Cameca,prefix='Plesovice',tst=c(337.13,0.18))
#' fit <- calibration(stand,oxide='UO2',omit=2)
#' calplot(stand,fit)
#' @export
calibration <- function(stand,oxide='UO2',omit=NULL){
    dat <- stand
    if (!is.null(omit)) dat$x <- stand$x[-omit]
    fit <- stats::optim(c(A=0,B=1),AB_misfit,stand=dat,oxide=oxide)
    hess <- stats::optimHess(fit$par,AB_misfit,stand=dat,oxide=oxide)
    out <- list()
    out$omit <- omit
    out$AB <- fit$par
    out$cov <- solve(hess)
    out$oxide <- oxide
    out$c64 <- stand$c64
    out$PbU <- stand$PbU
    out
}

AB_misfit <- function(AB,stand,oxide='UO2'){
    out <- 0
    snames <- names(stand$x)
    for (sname in snames){
        spot <- stand$x[[sname]]
        out <- out - LL_AB(AB,spot,oxide=oxide,c64=stand$c64)
    }
    out
}

LL_AB <- function(AB,spot,oxide='UO2',c64=18.7){
    p <- pars(spot,oxide=oxide)
    g <- get_gamma(B=AB['B'],p=p)
    b4 <- getPbLogRatio(p=p,g=g$Pb,num='4')
    bO <- log(p$cO/p$cU) + g$O*(p$tU-p$tO)
    a6 <- get_alpha(p=p,g=g$Pb,den='6')
    aU <- get_alpha(p=p,g=g$U,den='U')
    aO <- get_alpha(p=p,g=g$O,den='O')
    log_n6Ui <- AB[1] + AB[2]*(bO+aO-aU) + aU - a6 -
        log(1-exp(b4)*c64) - log(p$dU/p$d6) - g$Pb*(p$tU-p$t6)
    LL <- p$n6*log_n6Ui - (p$n6+p$nU)*log(1+exp(log_n6Ui))
    sum(LL)
}
