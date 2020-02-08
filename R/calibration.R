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
    bg <- get_bg(dat,oxide='UO2')
    fit <- stats::optim(c(A=0,B=1),AB_misfit,stand=stand,oxide=oxide,bg=bg)
    hess <- stats::optimHess(fit$par,AB_misfit,stand=stand,oxide=oxide,bg=bg)
    out <- list()
    out$omit <- omit
    out$AB <- fit$par
    out$cov <- solve(hess)
    out$oxide <- oxide
    out$c64 <- stand$c64
    out$PbU <- stand$PbU
    out
}

AB_misfit <- function(AB,stand=stand,oxide='UO2',bg=bg){
    snames <- names(stand$x)
    LL <- rep(0,length(snames))
    names(LL) <- snames
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        X <- log(p$O$c/p$U238$c) + bg$O[sname,'g']*(p$U238$t-p$O$t) +
            blank_correction(bg=bg$O[sname,c('b','g')],
                             bb=bg$blank[sname,'b'],
                             tt=p$O$t) -
            blank_correction(bg=bg$U[sname,c('b','g')],
                             bb=bg$blank[sname,'b'],
                             tt=p$U$t)
        RHS <- AB[1] + AB[2] * X
    }
    -sum(LL)
}

