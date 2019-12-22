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
    a <- get_alpha(AB=AB,p=p,g=g,c64=c64)
    # 1. get count ratios
    n6U <- exp(a$a6 + g$Pb*(p$t6-p$tU))*p$d6/p$dU
    # 2. get proportions
    nt <- length(p$tU)
    den <- n6U + 1
    theta <- cbind(n6U,1)/matrix(rep(den,2),ncol=2)
    colnames(theta) <- c('Pb206','U238')
    counts <- spot$counts[,c('Pb206','U238')]
    # 3. compute multinomial log-likelihood
    sum(counts*log(theta))
}
