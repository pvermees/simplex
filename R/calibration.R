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
    BG <- get_BG(dat,oxide='UO2')
    fit <- stats::optim(c(A=0,B=1),AB_misfit,stand=stand,oxide=oxide,BG=BG)
    hess <- stats::optimHess(fit$par,AB_misfit,stand=stand,oxide=oxide,BG=BG)
    out <- list()
    out$omit <- omit
    out$AB <- fit$par
    out$cov <- solve(hess)
    out$oxide <- oxide
    out$c64 <- stand$c64
    out$PbU <- stand$PbU
    out
}

AB_misfit <- function(AB,stand,oxide='UO2',BG){
    snames <- names(stand$x)
    out <- 0
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        cc <- get_cal_components(p=p,bg=BG[[sname]])
        # predict Pb206/U238 cps logratio (bp6U):
        X <- cc$bmOU - cc$dcOU + cc$bcO - cc$bcU
        RHS <- AB[1] + AB[2] * X
        b4corr <- log(1 - exp(cc$bdc46)*stand$c64)
        bp6U <- RHS + cc$dc6U - cc$bc6 + cc$bcU - b4corr
        # calculate the log-likelihood
        b6Ucounts <- bp6U+log(p$Pb206$d)-log(p$U238$d)
        n6 <- p$Pb206$n
        nU <- p$U238$n
        LL <- n6*b6Ucounts - (n6+nU)*log(1+exp(b6Ucounts))
        out <- out - sum(LL)
    }
    out
}

# bg = matrix of beta-gamma values for a single spot
get_cal_components <- function(p,bg){
    out <- list()
    tU <- p$U238$t
    # bm = measured beta (logratio of cps):
    out$bmOU <- log(p$O$c) - log(p$U238$c)
    out$bm46 <- bg['Pb204','b'] - bg['Pb206','b']
    # bc = blank correction:
    out$bcO <- blank_correct(bg=bg,tt=tU,mass='O')
    out$bcU <- blank_correct(bg=bg,tt=tU,mass='U238')
    out$bc4 <- blank_correct(bg=bg,tt=tU,mass='Pb204')
    out$bc6 <- blank_correct(bg=bg,tt=tU,mass='Pb206')
    # dc = drift correction
    out$dcOU <- bg['O','g']*(p$O$t-p$U238$t)
    out$dc46 <- bg['Pb206','g']*(p$Pb204$t-p$Pb206$t)
    out$dc6U <- bg['Pb206','g']*(p$Pb206$t-p$U238$t)
    # blank and drift corrected Pb204/Pb206 logratio
    out$bdc46 <- out$bm46 + out$dc46 + out$bc4 - out$bc6
    out
}
