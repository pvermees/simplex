#' @title create a calibration curve
#' @description fit a straight line to log[Pb/U] vs. log[UOx/U] data
#' @details Regresses a line through
#'     X=\eqn{^{238}}U\eqn{^{16}}O\eqn{_x}/\eqn{^{238}}U and
#'     Y=(\eqn{^{206}}Pb-\eqn{^{204}}Pb[\eqn{^{206}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{238}}U
#'
#' or through
#' 
#'     X=\eqn{^{232}}Th\eqn{^{16}}O\eqn{_x}/\eqn{^{232}}Th and
#'     Y=(\eqn{^{208}}Pb-\eqn{^{204}}Pb[\eqn{^{208}}Pb/
#'     \eqn{^{204}}Pb]\eqn{_{\circ}})/\eqn{^{232}}Th
#' @param stand a dataset of class \code{simplex}
#' @param oxide either \code{'UO'}, \code{'UO2'}, \code{'ThO'}, or
#'     \code{'ThO2'}
#' @param parent either \code{'U238'} or \code{'Th232'}
#' @param daughter either \code{'Pb206'} or \code{'Pb208'}
#' @param cD4 common \eqn{^{206}}Pb/\eqn{^{204}}Pb ratio (if
#'     \code{parent = 'U238'}) or common \eqn{^{208}}Pb/\eqn{^{204}}Pb
#'     ratio (if \code{parent = 'Th232'})
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
calibration <- function(stand,oxide='UO2',parent='U238',
                        daughter='Pb206',cD4=0,omit=NULL){
    dat <- stand
    if (!is.null(omit)) dat$x <- stand$x[-omit]
    B0G <- get_B0G(dat$x,parent=parent,oxide=oxide)
    init <- initAB(stand=stand,parent=parent,daughter=daughter,oxide=oxide)
    fit <- stats::optim(init,misfit_AB,stand=stand,oxide=oxide,
                        parent=parent,daughter=daughter,cD4=cD4,B0G=B0G)
    hess <- stats::optimHess(fit$par,misfit_AB,stand=stand,oxide=oxide,
                             parent=parent,daughter=daughter,cD4=cD4,B0G=B0G)
    out <- list()
    class(out) <- "calibration"
    out$parent <- parent
    out$daughter <- daughter
    out$oxide <- oxide
    out$cD4 <- cD4
    out$omit <- omit
    out$AB <- fit$par
    out$AB.cov <- solve(hess)
    out$DP <- stand$DP
    out$DP.cov <- stand$DP.cov
    out
}

initAB <- function(stand,parent='U238',daughter='Pb206',oxide='UO2'){
    X <- NULL
    Y <- NULL
    snames <- names(stand$x)
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],parent=parent,
                  daughter=daughter,oxide=oxide)
        bOP <- log(p[[p$oxide]]$c) - log(p[[p$parent]]$c)
        bDP <- log(p[[p$daughter]]$c) - log(p[[p$parent]]$c)
        X <- c(X,bOP)
        Y <- c(Y,bDP)
    }
    out <- lm(Y ~ X)$coef
    names(out) <- c('A','B')
    out
}

misfit_AB <- function(AB,stand,oxide='UO2',parent='U238',
                      daughter='Pb206',cD4=0,B0G){
    snames <- names(stand$x)
    out <- 0
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],parent=parent,
                  daughter=daughter,oxide=oxide)
        b0g <- B0G[[sname]]
        XY <- getCalXY(p=p,b0g=b0g,cD4=cD4)
        RHS <- AB[1] + AB[2] * XY$X
        bDPmc <- log(p[[p$daughter]]$c) - log(p[[p$parent]]$c)
        bDPpc <- RHS - XY$Y + bDPmc
        bDPpn <- bDPpc + log(p[[p$daughter]]$d) - log(p[[p$parent]]$d)
        LL <- LLbinom(bn=bDPpn,nnum=p[[p$daughter]]$n,nden=p[[p$parent]]$n)
        out <- out - sum(LL)
    }
    out
}
getCalXY <- function(p,b0g,cD4=0){
    out <- list()
    # get X
    bOPmc <- log(p[[p$oxide]]$c) - log(p[[p$parent]]$c)
    bdcorrOP <- A2Corr(p=p,b0g=b0g,num=p$oxide,den=p$parent)
    out$X <- bOPmc - bdcorrOP
    # get Y
    bDPmc <- log(p[[p$daughter]]$c) - log(p[[p$parent]]$c)
    bdcorr4D <- A2Corr(p=p,b0g=b0g,num='Pb204',den=p$daughter)
    b4D <- b0g['Pb204','b0'] - b0g[p$daughter,'b0'] - bdcorr4D
    b4corr <- log(1 - exp(b4D)*cD4)
    bdcorrDP <- A2Corr(p=p,b0g=b0g,num=p$daughter,den=p$parent)
    out$Y <- bDPmc - bdcorrDP + b4corr
    out
}
