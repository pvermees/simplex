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
    B0G <- get_B0G(dat,oxide=oxide)
    init <- initAB(stand=stand,oxide=oxide,parent=parent,daughter=daughter)
    fit <- stats::optim(init,misfit_AB,stand=stand,oxide=oxide,
                        parent=parent,daughter=daughter,cD4=cD4,B0G=B0G)
    hess <- stats::optimHess(fit$par,misfit_AB,stand=stand,oxide=oxide,
                             parent=parent,daughter=daughter,cD4=cD4,B0G=B0G)
    out <- list()
    out$parent <- parent
    out$daughter <- daughter
    out$cD4 <- cD4
    out$omit <- omit
    out$AB <- fit$par
    out$cov <- solve(hess)
    out$oxide <- oxide
    out$PbU <- stand$PbU
    out
}

initAB <- function(stand,oxide='UO2',parent='U238',daughter='Pb206'){
    X <- NULL
    Y <- NULL
    snames <- names(stand$x)
    for (sname in snames){
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        bOP <- log(p$O$c) - log(p[[parent]]$c)
        bDP <- log(p[[daughter]]$c) - log(p[[parent]]$c)
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
        p <- pars(spot=stand$x[[sname]],oxide=oxide)
        b0g <- B0G[[sname]]
        out <- out + misfit_A(A=AB[1],B=AB[2],p=p,b0g=b0g,
                              parent=parent,daughter=daughter,cD4=cD4)
    }
    out
}
misfit_A <- function(A,B,p,b0g,parent='U238',daughter='Pb206',
                     cD4=0,deriv=FALSE){
    XY <- getCalXY(p=p,b0g=b0g,parent=parent,daughter=daughter,cD4=cD4)
    RHS <- A + B * XY$X
    bDPmc <- log(p[[daughter]]$c) - log(p[[parent]]$c)
    bDPpc <- RHS - XY$Y + bDPmc
    bDPpn <- bDPpc + log(p[[daughter]]$d) - log(p[[parent]]$d)
    LL <- LLbinom(bn=bDPpn,nnum=p[[daughter]]$n,nden=p[[parent]]$n)
    if (deriv)
        out <- dAdB(X=XY$X,bDPpc=bDPpc,nD=p[[daughter]]$n,nP=p[[parent]]$n)
    else
        out <- (-sum(LL))
    out
}
getCalXY <- function(p,b0g,parent='U238',daughter='Pb206',cD4=0){
    out <- list()
    # get X
    bOPmc <- log(p$O$c) - log(p[[parent]]$c)
    bdcorrPO <- A2Corr(p=p,b0g=b0g,num='O',den=parent)
    out$X <- bOPmc - bdcorrPO
    # get Y
    bDPmc <- log(p[[daughter]]$c) - log(p[[parent]]$c)
    bdcorrDP <- A2Corr(p=p,b0g=b0g,num=daughter,den=parent)
    bdcorr4D <- A2Corr(p=p,b0g=b0g,num='Pb204',den=daughter)
    b4D <- b0g['Pb204','b0'] - b0g[daughter,'b0'] - bdcorr4D
    b4corr <- log(1 - exp(b4D)*cD4)
    out$Y <- bDPmc - bdcorrDP + b4corr
    out
}
dAdB <- function(X,bDPpc,nD,nP){
    dbDPpcdA <- 1
    dbDPpcdB <- X
    num <- exp(bDPpc)
    den <- 1 + exp(bDPpc)
    dnumdA <- exp(bDPpc)*dbDPpcdA
    dnumdB <- exp(bDPpc)*dbDPpcdB
    dLdA <- dbDPpcdA*(nD-(nD+nP)*num/den)
    d2LdA2 <- -dbDPpcdA*dnumdA*(nD+nP)/(den^2)
    d2LdAdB <- -dbDPpcdA*dnumdB*(nD+nP)/(den^2)
    -sum(d2LdAdB)/sum(d2LdA2)
}
