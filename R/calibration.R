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
    out <- 0
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        tU <- p$U238$t
        # bm = measured beta (logratio of cps):
        bmOU <- log(p$O$c) - log(p$U238$c)
        bm46 <- bg$Pb204[sname,'b'] - bg$Pb206[sname,'b']
        # bc = blank correction:
        bcO <- blank_correct(bg=bg,sname=sname,tt=tU,mass='O')
        bcU <- blank_correct(bg=bg,sname=sname,tt=tU,mass='U238')
        bc4 <- blank_correct(bg=bg,sname=sname,tt=tU,mass='Pb204')
        bc6 <- blank_correct(bg=bg,sname=sname,tt=tU,mass='Pb206')
        # dc = drift correction
        dcOU <- bg$O[sname,'g']*(p$O$t-p$U238$t)
        dc46 <- bg$Pb206[sname,'g']*(p$Pb204$t-p$Pb206$t)
        dc6U <- bg$Pb206[sname,'g']*(p$Pb206$t-p$U238$t)
        # infer bm6U:
        X <- bmOU - dcOU + bcO - bcU
        RHS <- AB[1] + AB[2] * X
        b4corr <- log(1 - exp(bm46 + dc46 + bc4 - bc6)*stand$c64)
        bm6U <- RHS + dc6U - bc6 + bcU - b4corr
        # calculate the log-likelihood
        b6Ucounts <- bm6U+log(p$Pb206$d)-log(p$U238$d)
        n6 <- p$Pb206$n
        nU <- p$U238$n
        LL <- n6*b6Ucounts - (n6+nU)*log(1+exp(b6Ucounts))
        out <- out - sum(LL)
    }
    out
}

